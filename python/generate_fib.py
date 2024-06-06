MAX_LEN = 1 << 27

f1 = "b"
f2 = "a"

while len(f1) < MAX_LEN:
    f1, f2 = f2, f2 + f1

print(f1[:MAX_LEN])

