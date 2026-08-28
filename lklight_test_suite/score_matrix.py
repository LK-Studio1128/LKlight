import subprocess, time, csv
ROOT = "/Users/luoxiaowen/Desktop/LKDock/byi/LKlight"
BIN = f"{ROOT}/target/release/LKlight"
OUT = "/tmp/lktest"
methods = "dfire fastdfire dfire2 dna ddna mj3h pydock cpydock sd vdw pisa sipper tobi".split()
cases = [
    ("1azp", f"{ROOT}/tests/1azp/1azp_receptor.pdb", f"{ROOT}/tests/1azp/1azp_ligand.pdb"),
    ("2oob", f"{ROOT}/tests/2oob/2oob_receptor.pdb", f"{ROOT}/tests/2oob/2oob_ligand.pdb"),
    # 场景案例也纳入 score 矩阵（代表性函数集）
    ("p53DNA", "/tmp/lktest/cases/1DIZ_p53.pdb", "/tmp/lktest/cases/1DIZ_DNA.pdb"),
    ("AbLyso", "/tmp/lktest/cases/1VFB_antibody.pdb", "/tmp/lktest/cases/1VFB_lysozyme.pdb"),
    ("AbHIVpep", "/tmp/lktest/cases/1DQJ_antibody.pdb", "/tmp/lktest/cases/1DQJ_peptide.pdb"),
    ("RBD_ACE2", "/tmp/lktest/cases/6M0J_ACE2.pdb", "/tmp/lktest/cases/6M0J_RBD.pdb"),
]
rows = []
for name, rec, lig in cases:
    for m in methods:
        t0 = time.time()
        try:
            p = subprocess.run([BIN, "score", rec, lig, m], capture_output=True, text=True, timeout=120)
            score = "ERR"
            for line in p.stdout.splitlines():
                if "Score" in line and ":" in line:
                    try:
                        score = float(line.split(":")[1])
                    except ValueError:
                        pass
        except subprocess.TimeoutExpired:
            score = "TIMEOUT"
        dt = (time.time() - t0) * 1000
        rows.append([name, m, score, f"{dt:.1f}"])
        print(rows[-1], flush=True)
with open(f"{OUT}/scores/matrix.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["case", "method", "score", "time_ms"])
    w.writerows(rows)
print("saved", f"{OUT}/scores/matrix.csv")
