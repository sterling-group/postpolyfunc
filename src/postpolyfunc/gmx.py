from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple, Union
import shutil
import gmxapi as gmx


@dataclass
class StepOutput:
    name: str
    files: Dict[str, str] = field(default_factory=dict)  # logical key -> absolute path
    meta: Dict[str, Any] = field(default_factory=dict)   # anything else (timings, logs, params)


class GmxAPI:
    """
    Thin wrapper around gmxapi + an orchestrator to run multiple steps in sequence.
    You keep filenames internal; only parameters like box size come from CLI.
    """

    def __init__(
        self,
        executable: str = "gmx_mpi",
        workdir: Union[str, Path, None] = None,
        args: Any = None,  # argparse.Namespace with things like .box, if you have it
    ):
        self.executable = executable
        self.workdir = Path(workdir) if workdir else Path.cwd()
        self.args = args
        self.workdir.mkdir(parents=True, exist_ok=True)

        # Registry maps step name -> bound method
        self._registry: Dict[str, Callable[..., StepOutput]] = {
            "create_box": self.create_boxed_structure_step,
            # Add more as you implement them:
            "grompp_em": self.grompp_em_step,   # example (energy-min input)
            "mdrun_em": self.mdrun_em_step,     # example (energy-min MD)
            # "solvate": self.solvate_step,
            # "ions": self.genion_step,
            # "nvt": self.nvt_step,
            # ...
        }

    # ---------- Low-level helper ----------

    def _run_cmd(
        self,
        arguments: List[str],
        input_files: Optional[Dict[str, str]] = None,
        output_files: Optional[Dict[str, str]] = None,
    ) -> Dict[str, str]:
        """
        Thin wrapper around gmx.commandline_operation using keyword args (prevents signature errors).
        Returns a dict of realized output file paths keyed by their CLI flag (e.g., "-o").
        """
        op = gmx.commandline_operation(
            command=[self.executable],
            arguments=arguments,
            input_files=input_files or {},
            output_files=output_files or {},
        )
        op.run()

        realized: Dict[str, str] = {}
        for flag in (output_files or {}):
            realized[flag] = str(Path(op.output.file[flag].result()).resolve())
        return realized

    # ---------- Steps (return StepOutput) ----------

    def create_boxed_structure_step(
        self,
        input_pdb: str,
        box: Optional[Iterable[float]] = None,
        boxtype: str = "cubic",
        center: bool = False,  # upstream already centers; keep False by default
        outname: str = "solute_boxed.gro",
    ) -> StepOutput:
        """
        Solute-only box creation (editconf). Output name is internal/stable for downstream steps.
        """
        if box is None and self.args is not None:
            box = getattr(self.args, "box", None)
        if box is None:
            box = (12.0, 12.0, 12.0)
        box = list(box)
        if len(box) != 3:
            raise ValueError(f"box must be 3 numbers, got {box}")

        out_path = (self.workdir / outname).resolve()
        arguments = ["editconf"]
        if center:
            arguments.append("-c")
        arguments += ["-bt", boxtype, "-box", str(box[0]), str(box[1]), str(box[2])]

        realized = self._run_cmd(
            arguments=arguments,
            input_files={"-f": str(input_pdb)},
            output_files={"-o": str(out_path)},
        )

        produced = Path(realized["-o"]).resolve()
        if produced != out_path:
            out_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(produced), str(out_path))

        return StepOutput(
            name="create_box",
            files={"gro": str(out_path)},
            meta={"boxtype": boxtype, "box": tuple(float(x) for x in box), "center": center},
        )

    def grompp_em_step(
        self,
        gro: str,
        top: str,
        mdp: str,
        out_tpr: str = "em.tpr",
    ) -> StepOutput:
        """
        Example: pre-process for energy minimization.
        Requires: structure (.gro), topology (.top/.itp included), and an EM .mdp.
        """
        out_tpr_path = (self.workdir / out_tpr).resolve()
        realized = self._run_cmd(
            arguments=["grompp"],
            input_files={
                "-f": str(mdp),
                "-c": str(gro),
                "-p": str(top),
            },
            output_files={"-o": str(out_tpr_path)},
        )
        return StepOutput(
            name="grompp_em",
            files={"tpr": realized["-o"]},
            meta={},
        )

    def mdrun_em_step(
        self,
        tpr: str,
        deffnm: str = "em",
    ) -> StepOutput:
        """
        Example: run energy minimization (mdrun).
        Produces em.gro, em.edr, em.log, em.trr (names hidden behind deffnm).
        """
        # With gmxapi, you can either pass explicit outputs or rely on -deffnm.
        # Here we use -deffnm to keep the file family consistent.
        realized = self._run_cmd(
            arguments=["mdrun", "-deffnm", deffnm],
            input_files={"-s": str(tpr)},
            output_files={
                # Ask gmxapi to track primary outputs you care about:
                "-c": str((self.workdir / f"{deffnm}.gro").resolve()),
                "-e": str((self.workdir / f"{deffnm}.edr").resolve()),
                "-g": str((self.workdir / f"{deffnm}.log").resolve()),
                "-o": str((self.workdir / f"{deffnm}.trr").resolve()),
            },
        )
        # Normalize keys for consumers
        return StepOutput(
            name="mdrun_em",
            files={
                "gro": realized.get("-c"),
                "edr": realized.get("-e"),
                "log": realized.get("-g"),
                "trr": realized.get("-o"),
            },
            meta={"deffnm": deffnm},
        )

    # ---------- Orchestrator ----------

    def orchestrate(
        self,
        steps: Iterable[str],
        *,
        inputs: Dict[str, Any],
        overrides: Optional[Dict[str, Dict[str, Any]]] = None,
        strict: bool = True,
    ) -> Dict[str, StepOutput]:
        """
        Run a sequence of registered steps, passing outputs forward.

        Parameters
        ----------
        steps : list of step names in execution order.
        inputs : dict of initial inputs (e.g., {"input_pdb": "...", "top": "...", "mdp": "..."})
        overrides : optional per-step kwargs to override defaults
                    e.g., { "create_box": {"box": (10,10,12)}, "mdrun_em": {"deffnm": "min"} }
        strict : if True, raise on unknown step; else skip with warning.

        Returns
        -------
        dict : step name -> StepOutput
        """
        ctx: Dict[str, Any] = dict(inputs)  # working context (files + params)
        results: Dict[str, StepOutput] = {}
        overrides = overrides or {}

        for name in steps:
            fn = self._registry.get(name)
            if fn is None:
                if strict:
                    raise KeyError(f"Unknown step: {name}")
                else:
                    print(f"[WARN] Skipping unknown step: {name}")
                    continue

            # Build call kwargs from ctx + specific overrides
            kw = {}
            # Naive but effective: pass only kwargs that the step likely uses
            # You can keep it simple and rely on Python to raise if missing.
            kw.update(ctx)
            kw.update(overrides.get(name, {}))

            # Run
            out = fn(**kw)  # type: ignore[arg-type]
            results[name] = out

            # Promote produced files into ctx using predictable keys
            # so downstream steps can reference them easily.
            # Example conventions:
            if name == "create_box":
                ctx["gro"] = out.files["gro"]
            elif name == "grompp_em":
                ctx["tpr"] = out.files["tpr"]
            elif name == "mdrun_em":
                # new minimized structure:
                if out.files.get("gro"):
                    ctx["gro"] = out.files["gro"]

            # You can also inject `meta` if downstream needs them
            ctx.update({f"{name}__meta": out.meta})

        return results
