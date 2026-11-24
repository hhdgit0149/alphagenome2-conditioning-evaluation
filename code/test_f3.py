import sys, importlib.util
spec = importlib.util.spec_from_file_location("figs", r"D:/ar_pilot/69_figures.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
cases = [(0.6154985596768819, "0.615"), (0.2124651037403786, "0.212"),
         (0.7160, "0.716"), (-0.024690347395356384, "-0.025"),
         (0.6070702430968664, "0.607")]
bad = [(x, m.f3(x), e) for x, e in cases if m.f3(x) != e]
print("cases", len(cases), "failures", bad)
assert not bad
print("f3 OK")
