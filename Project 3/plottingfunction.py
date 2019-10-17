import numpy as np
import matplotlib.pyplot as plt


#plotting the function with different limit values to see which limit to choose,
# we know that e^(-2*lam) = 0 for a good estimate lam  
def function(lam):
    return np.exp(-2*lam)

lam = np.linspace(1, 10, 10)
print(lam)

y = function(lam)

plt.plot(y)
plt.title('Checking the limit approximation')
plt.xlabel('lambda')
plt.ylabel('exp(-2*lambda)')
plt.show()