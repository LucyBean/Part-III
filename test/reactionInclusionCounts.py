import json
import matplotlib.pyplot as plt

fp = "FACOAL161t2pp EFMs.json"

with open(fp, "r") as f:
    s = f.read()
    data = json.loads(s)
    
# Convert to sets of reactions
reacSets = [set(k.keys()) for k in data]

# Count frequency of each reaction
reacFreq = {}
for rs in reacSets:
    for r in rs:
        if r not in reacFreq:
            reacFreq[r] = 0
        reacFreq[r] += 1
        
freqs = sorted(reacFreq.values())

plt.plot(freqs)
plt.ylabel("Frequency of inclusion")
plt.xlabel("Reaction")
plt.show()
