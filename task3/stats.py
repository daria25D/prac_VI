import statistics as stat
import matplotlib.pyplot as plt
import seaborn as sns
import sys

with open(sys.argv[1]) as f:
        numbers = []
        for line in f:
            numbers.append(float(line))
sns.set(style='white', palette='muted', color_codes=True)
plt.xlabel("1-F")
plt.ylabel("Occurences")
sns.distplot(numbers, bins=20, kde=False, color="b")
plt.tight_layout()
plt.show()