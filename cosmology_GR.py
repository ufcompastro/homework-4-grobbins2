#Homework 4 worked on by Grady Robbins
from gaussxw import gaussxw,gaussxwab
import matplotlib.pyplot as plt
import numpy as np

def H(z):
    return (QM*((1+z)**3) +QR*((1+z)**2) + QA)**.5
def r(z):
    return C / H(z) #defining functions for integration
z1 = 0
C = 3000.
QM = 0.3
QA = 0.7
QR = 0
z2f = 10 # defining variables
print("Through linear integration from z1=0 to z2=10, Angular Diameter is:")
Int_lin = 0
Nlin = 50
DA_lin = np.zeros((11),float)

for z2 in range(z1, z2f+1): # using for loop to create an array of linear integrations
    Int_lin = 0
    for z in np.linspace(z1,z2,Nlin):
        h = (z2 - z1)/Nlin
        A_rect = h*r(z)
        Int_lin += A_rect
    DA_lin[z2] = Int_lin/(1+z2)
print(DA_lin)

print("Through trapezoidal integration from z1=0 to z2=10, Angular Diameter is:")
Ntrap = 50
DA_trap = np.zeros((11),float) 
for z2 in range(z1, z2f+1): # repeating for trapezoidal
    int_trap = 0
    for k in range(Ntrap):
        h = (z2 - z1)/Ntrap
        int_trap += r(z1+k*h)
    DA_trap[z2] = (h*((.5*r(z1) + .5*r(z2) + int_trap)))/(1+z2)
print(DA_trap)

print("Through simpson's rule from z1=0 to z2=10, Angular Diameter is:")
DA_simpson = np.zeros((11),float)
Nsimp = 50
for z2 in range(z1, z2f+1): # repeating for simpson's rule
    Int_even = 0
    Int_odd = 0
    for k in range(1,Nsimp,2):
        h = (z2 - z1)/Nsimp
        Int_odd += r(z1+k*h)
    for k in range(2,Nsimp,2):
        h = (z2 - z1)/Nsimp
        Int_even += r(z1+k*h)
    Int_simpson = (h/3)*(r(z1) + r(z2) + 4*Int_odd + 2*Int_even)
    DA_simpson[z2] = Int_simpson/(1+z2)
print(DA_simpson)

print("Through Gaussian Quadrature from z1=0 to z2=10, Angular Diameter is:")
Ngauss = 50
DA_gauss = np.zeros((11),float)
for z2 in range(z1,z2f+1): # repeating for gaussian quadrature
    a = z1
    b = z2
    x,w = gaussxw(Ngauss)
    xprime,wprime = 0.5*(b-a)*x + 0.5*(b+a) , 0.5 * (b-a)*w
    Int_gauss = 0
    for k in range(Ngauss):
        Int_gauss = Int_gauss + wprime[k]*r(xprime[k])
    DA_gauss[z2] = Int_gauss/(1+z2)
print(DA_gauss)
print("Through np.trapz from z1=0 to z2=10, Angular Diameter is:") # creating np.trapz specific arrays for different integration technique
nptrapz = np.zeros((11),float)
z_shift = np.linspace(0,11,11)
ang_z = np.zeros((11),float)
for n in range(0,11):
    ang_z[n] = r(n) 
for k in range(0,11):
    nptrapz[k] = np.trapz(ang_z[0:k+1],z_shift[0:k+1]) / (1 + k)
print(nptrapz)
plt.plot(z_shift,nptrapz,color='green', marker='o', markersize=8, linestyle=':', label = 'np.trapz') # adding graphs with unique colors and variables 
plt.plot(z_shift,DA_lin,color='blue', marker='o', markersize=8, linestyle=':', label = 'linear')
plt.plot(z_shift,DA_trap,color='orange', marker='o',markersize=8, linestyle=':', label = 'trapezoidal')
plt.plot(z_shift,DA_simpson,color='red', marker='o',markersize=8, linestyle=':', label = 'simpson')
plt.plot(z_shift,DA_gauss,color='cyan', marker='o',markersize=5, linestyle=':', label = 'gaussian')
plt.xlabel('Redshift (z)')
plt.ylabel('Angular Diameter')
plt.legend() #creating legend, labels, and saving graph
plt.savefig('robbins_hw4.png', dpi=300)
