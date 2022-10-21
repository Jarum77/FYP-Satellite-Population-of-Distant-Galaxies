from mpmath import mp
import NDpredict as ndp
mp.dps = 50
print(mp.quad(lambda x: mp.exp(-x**2), [-mp.inf, mp.inf]) ** 2)