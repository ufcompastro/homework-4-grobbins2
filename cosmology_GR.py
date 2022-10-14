#Homework 4 worked on by Grady Robbins
from gaussxw import gaussxw,gaussxwab
import numpy as np

def H(z):
    return (QM*((1+z)**3) +QR*((1+z)**2) + QA)**.5
def r(z):
    return C / H(z)
z1 = 0
C = 3000.
QM = 0.3
QA = 0.7
QR = 0
z2f = 10
print("Through linear integration, Angular Diameter is:")
Int_lin = 0
Nlin = 10
DA_lin = np.zeros((1,11),float)

for z2 in range(z1, z2f+1):
    Int_lin = 0
    for z in np.linspace(z1,z2,Nlin):
        h = (z2 - z1)/Nlin
        A_rect = h*r(z)
        Int_lin += A_rect
    DA_lin[0,z2] = Int_lin/(1+z2)
print(DA_lin)

print("Through trapezoidal integration, Angular Diameter is:")
Ntrap = 10
DA_trap = np.zeros((1,11),float)
for z2 in range(z1, z2f+1):
    int_trap = 0
    for k in range(Ntrap):
        h = (z2 - z1)/Ntrap
        int_trap += r(z1+k*h)
    DA_trap[0,z2] = (h*((.5*r(z1) + .5*r(z2) + int_trap)))/(1+z2)
print(DA_trap)

print("Through simpson's rule, Angular Diameter is:")
DA_simpson = np.zeros((1,11),float)
Nsimp = 10
for z2 in range(z1, z2f+1):
    Int_even = 0
    Int_odd = 0
    for k in range(1,Nsimp,2):
        h = (z2 - z1)/Nsimp
        Int_odd += r(z1+k*h)
    for k in range(2,Nsimp,2):
        h = (z2 - z1)/Nsimp
        Int_even += r(z1+k*h)
    Int_simpson = (h/3)*(r(z1) + r(z2) + 4*Int_odd + 2*Int_even)
    DA_simpson[0,z2] = Int_simpson/(1+z2)
print(DA_simpson)

print("Through Gaussian Quadrature, Angular Diameter is:")
Ngauss = 1
DA_gauss = np.zeros((1,11),float)
#for z2 in range(z1,z2f+1):
a = 0
b = 3
x,w = gaussxwab(Ngauss,a,b)
#xprime,wprime = 0.5*(b-a)*x + 0.5*(b+a) , 0.5 * (b-a)*w
Int_gauss = 0
for k in range(Ngauss):
    Int_gauss = Int_gauss + w[k] + r(x[k])
print(Int_gauss/2)
#DA_gauss[0,z2] = Int_gauss/(1+z2f)
#print(DA_gauss)
