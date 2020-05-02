import matplotlib.pyplot as plt
import numpy as np
import math


plot1 = plt.figure(figsize=(2,6.2832),dpi=80)

plot1.add_subplot(2,2,2)
x = np.arange(0,np.pi*2,0.1)
y = np.sin(x)
plt.plot(x,y)
plt.legend(["sin"])
