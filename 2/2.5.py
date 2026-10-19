#19.10.
s = input()
t = input()

print(sum(a != b for a, b in zip(s, t)))
