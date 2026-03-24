def pytest_addoption(parser):
    parser.addoption(
        "--run-gpu-smoke",
        action="store_true",
        default=False,
        help="Run smoke tests that require GPU.",
    )
    parser.addoption(
        "--smoke-task",
        action="append",
        default=[],
        help="Run only named smoke tasks (can be repeated).",
    )
    parser.addoption(
        "--smoke-timeout-multiplier",
        action="store",
        type=float,
        default=1.0,
        help="Multiply per-case smoke timeout values by this factor.",
    )
