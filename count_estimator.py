"""
Each mutation may cause a difference between the original and mutated DNA string.

To estimate the total number of mutations, we can count the total number of differences.
"""

# s is the original DNA string
s = input()

# t is the mutated DNA string
t = input()

n = len(s)
m = len(t)
c = n if n < m else m

# count is the total number of bases that differ between the original and mutated DNA
# strings.
count = abs(n - m)
for i in range(c):
    if s[i] != t[i]:
        count += 1

print(count)
