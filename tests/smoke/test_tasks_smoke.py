import signal

import pytest

from tests.smoke.manifest import SMOKE_CASES
from tests.smoke.validators import run_validations


@pytest.mark.smoke
@pytest.mark.parametrize("case", SMOKE_CASES, ids=[case["name"] for case in SMOKE_CASES])
def test_task_smoke(case, tmp_path, request):
    pytest.importorskip("ribbon")
    task_mod = pytest.importorskip("ribbon_tasks.tasks")

    only_tasks = request.config.getoption("--smoke-task") or []
    if only_tasks and case["name"] not in only_tasks and case["class_name"] not in only_tasks:
        pytest.skip("Not selected by --smoke-task.")

    if case.get("requires_gpu", False) and not request.config.getoption("--run-gpu-smoke"):
        pytest.skip("GPU smoke disabled. Use --run-gpu-smoke to enable.")

    case_dir = case["case_dir"]
    expected_dir = case_dir / "expected"
    out_dir = tmp_path / "out"
    out_dir.mkdir(parents=True, exist_ok=True)

    task_class = getattr(task_mod, case["class_name"])
    kwargs = case["kwargs_builder"](case_dir, out_dir)
    task = task_class(**kwargs)

    timeout_multiplier = request.config.getoption("--smoke-timeout-multiplier")
    timeout_s = max(1, int(case["timeout_s"] * timeout_multiplier))

    def _timeout_handler(_signum, _frame):
        raise TimeoutError(f"Smoke test for {case['name']} timed out after {timeout_s}s")

    signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(timeout_s)
    try:
        task.run()
    finally:
        signal.alarm(0)

    run_validations(case["validations"], out_dir, expected_dir)
