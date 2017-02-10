import json
import matplotlib.pyplot as plt
import numpy as np

# Load data if don't already have it
if "data" not in locals():
    with open("FACOAL161t2pp counts.json","r") as f:
        s = f.read()
        data = json.loads(s)
        
# Extract fields
times = data["times"]
dcounts = data["dcounts"]
icounts = data["icounts"]
ucounts = data["ucounts"]

def getDerivs(counts, frame, smooth=True):
    derivs = []
    # These are done as the difference between current count and
    # count in 10 steps time over the time difference
    for i in range(1,len(counts)-frame):
        dc = counts[i+frame] - counts[i]
        dt = times[i+frame] - times[i]
        derivs.append(dc/dt)
        
    if smooth:
        smoothDerivs = [derivs[0]]
        smoothing = 0.2
        for i in range(1,len(derivs)):
            smoothDerivs.append(smoothDerivs[i-1] * (1-smoothing) + derivs[i] * smoothing)
        return np.array(smoothDerivs)
    else:
        return np.array(derivs)
    
def graphConvert(count):
    return np.log(getDerivs(count,frame)+1)/np.log(10)

# Convert to "derivatives"
frame = 100
    
x = times[1:-frame]
ip, = plt.plot(x, graphConvert(icounts), label="Infeasible", color='red')
dp, = plt.plot(x, graphConvert(dcounts), label="Duplicate", color='green')
up, = plt.plot(x, graphConvert(ucounts), label="Unique", color="blue")
plt.legend(handles=[ip, dp, up], bbox_to_anchor=(0.4, 1))
plt.xlabel("Measurement index")
#plt.ylabel("log(Discovery rate,10)")
plt.ylabel("Discovery rate")
plt.show()
