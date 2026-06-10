import subprocess
import sys
from datetime import datetime
from pathlib import Path


def run_script(script_name, log_file):
    log_path = Path(log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 80, flush=True)
    print(f"Start running: {script_name}", flush=True)
    print(f"Log file: {log_file}", flush=True)
    print(f"Time: {datetime.now()}", flush=True)
    print("=" * 80, flush=True)

    with open(log_path, "w") as f:
        f.write("=" * 80 + "\n")
        f.write(f"Start running: {script_name}\n")
        f.write(f"Time: {datetime.now()}\n")
        f.write("=" * 80 + "\n")
        f.flush()

        result = subprocess.run(
            [sys.executable, "-u", script_name],
            stdout=f,
            stderr=subprocess.STDOUT,
            text=True,
        )

        f.write("\n" + "=" * 80 + "\n")
        if result.returncode != 0:
            f.write(f"ERROR: {script_name} failed with return code {result.returncode}\n")
            f.write(f"Time: {datetime.now()}\n")
            f.write("Stop here. The next script will NOT run.\n")
        else:
            f.write(f"Finished successfully: {script_name}\n")
            f.write(f"Time: {datetime.now()}\n")
        f.write("=" * 80 + "\n")
        f.flush()

    if result.returncode != 0:
        print("=" * 80, flush=True)
        print(f"ERROR: {script_name} failed with return code {result.returncode}", flush=True)
        print(f"See log: {log_file}", flush=True)
        print(f"Time: {datetime.now()}", flush=True)
        print("Stop here. The next script will NOT run.", flush=True)
        print("=" * 80, flush=True)
        sys.exit(result.returncode)

    print("=" * 80, flush=True)
    print(f"Finished successfully: {script_name}", flush=True)
    print(f"See log: {log_file}", flush=True)
    print(f"Time: {datetime.now()}", flush=True)
    print("=" * 80, flush=True)


if __name__ == "__main__":
    run_script(
        "3_1_adata_X_sparse260518.py",
        "log/3_1_adata_X_sparse260518.log",
    )

    run_script(
        "3_2_adata_layers_sparse260518.py",
        "log/3_2_adata_layers_sparse260518.log",
    )

    print("=" * 80, flush=True)
    print("All scripts finished successfully.", flush=True)
    print(f"Time: {datetime.now()}", flush=True)
    print("=" * 80, flush=True)