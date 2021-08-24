import sys
import glob

def disassemble(lines):
    comments, count = -1, 0
    ret = []
    for l in lines:
        if len(l) > 1:
            if l[0] == '#':
                count += 1
            else:
                if comments == -1:
                    comments = count
                assert comments == count
                assert comments % 2 == 0
                row = list(map(float, l.strip().split()))
                assert len(row) == comments
                separated = [row[:10], row[10:]]
                ret.append(separated)
    return ret

def calc_diff(a, b):
    assert len(a) == len(b)
    max_abs_diff = 0
    max_rel_diff = 0
    for i in range(len(a)):
        abs_diff = abs(a[i] - b[i])
        if abs(a[i]) < 1e-15:
            rel_diff = 1e9
        else:
            rel_diff = abs((b[i] - a[i]) / a[i])
        if abs_diff > max_abs_diff:
            max_abs_diff = abs_diff
        if rel_diff > max_rel_diff:
            max_rel_diff = rel_diff
    return [max_abs_diff, max_rel_diff]

#a and b are file paths of outputs we want to compare
def check(a, b):
    assert a != b, "arguments must differ!"
    lines_a, lines_b = None, None
    try:
        with open(a) as f:
            lines_a = f.readlines()
    except IOError:
        raise Exception("First file not accessible")

    try:
        with open(b) as f:
            lines_b = f.readlines()
    except IOError:
        raise Exception("Second file not accessible")

    d_a = disassemble(lines_a)
    d_b = disassemble(lines_b)

    max_abs_diff, max_rel_diff = 0, 0
    assert len(d_a) == len(d_b)
    for i in range(len(d_a)):
        assert d_a[i][0] == d_b[i][0]
        abs_diff, rel_diff = calc_diff(d_a[i][1], d_b[i][1])
        if abs_diff > max_abs_diff:
            max_abs_diff = abs_diff
        if rel_diff > max_rel_diff:
            max_rel_diff = rel_diff
    print("number of observations:", len(d_a))
    print("number of channels:", int(len(d_a[0][1]) / 2))
    print("max abs diff: {}, max rel dif: {}".format(max_abs_diff,
    max_rel_diff))

if __name__ == "__main__":
    with open('aux/submission_index') as f:
        lines = f.readlines()
        assert len(lines) == 1
        index = int(lines[0])
    print("Index of the last submission: {}\n".format(index))
    files = glob.glob("*/submissions/rad-{}.tab".format(index))
    files.sort()
    for f in files:
        num = f[0:3]
        print("Test number: {}".format(num))
        check("{}/rad-785-798.tab".format(num),
        "{}/submissions/rad-{}.tab".format(num, index)) 
        print("-------------------------------------\n")