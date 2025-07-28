from pathlib import Path
import pytest
import shutil
import subprocess
import sys
import tempfile


PREFIX = Path(__file__).parent.parent


def run(path, args=None, replace=None):
    """Run example script."""

    if args is None:
        args = ""

    path = PREFIX / f"examples/{path}"
    if replace is not None:
        tmp = tempfile.TemporaryDirectory()
        shutil.copytree(path.parent, tmp.name, dirs_exist_ok=True)

        path = Path(tmp.name) / path.name
        with path.open() as f:
            content = f.read()

        for k, v in replace.items():
            content = content.replace(k, v)

        with path.open("w") as f:
            f.write(content)

    command = f"{sys.executable} {path} {args}"
    r = subprocess.run(command, shell=True, capture_output=True)
    if r.returncode != 0:
        print(r.stdout.decode())
        raise RuntimeError(r.stderr.decode())

@pytest.mark.example
def test_benchmark_backward():
    """Test the benchmark/backward example."""

    replace = { "N = 1000000": "N = 1000" }
    run("benchmark/backward.py", replace=replace)

@pytest.mark.example
def test_benchmark_forward():
    """Test the benchmark/forward example."""

    replace = { "N = 1000000": "N = 1000" }
    run("benchmark/forward.py", replace=replace)

@pytest.mark.example
@pytest.mark.requires_calzone
def test_mixed():
    """Test the mixed example."""

    run("mixed/run-mixed.py -n 1000")

@pytest.mark.example
@pytest.mark.requires_calzone
def test_mixed_forward():
    """Test the mixed/forward example."""

    run("mixed/run-forward.py", args="-n 1000")

@pytest.mark.example
def test_transport_backward():
    """Test the transport/backward example."""

    replace = { "N = 1000000": "N = 1000" }
    run("transport/backward.py", replace=replace)

@pytest.mark.example
def test_transport_forward():
    """Test the transport/forward example."""

    replace = { "N = 1000000": "N = 1000" }
    run("transport/forward.py", replace=replace)

@pytest.mark.example
@pytest.mark.requires_matplotlib
def test_processes():
    """Test the processes example."""

    replace = { "plot.show()": "" }
    run("processes.py", replace=replace)
