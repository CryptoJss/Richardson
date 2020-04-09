import math as m
import sympy  as s
import matplotlib.pyplot as plt
from scipy import integrate

class Help:
    def graph(self,x1,fx1,title,x2,fx2):
        plt.plot(x1, fx1, label='F(x1)')
        plt.plot(x2, fx2, label='F(x2)')
        plt.title(title)
        plt.ylabel('f(x)')
        plt.xlabel('Xn')
        plt.legend()
        plt.show()

    def graphSimpson(self,x,fx,title):
        plt.plot(x, fx, label='F(x)')
        plt.title(title)
        plt.ylabel('f(x)')
        plt.xlabel('Xn')
        plt.legend()
        plt.show()

    def TrapezoidWidth(self,a,b,n):
        H = (b-a)/n
        return round(H,6)

    def Error(self,Ireal,Iaproxima):
        e = ((Ireal - Iaproxima)/Ireal)*100
        return round(abs(e),4)

class Trapezoidal(Help):
    def RichardsonTrape(self,fx,n,x):
        sfx = 0
        for i in range(1,n):
            sfx+=fx[i]
        I = (x[n]-x[0])*(fx[0]+2*sfx+fx[n])/(2*n)
        return round(I,6)

    def balance(self,I1,I2):
        Iaprox = 4/3*(I2) - 1/3*(I1)
        return round(Iaprox,6)

    def body(self):
        x = s.Symbol('x')
        n1 = int(input("Con Cuantos terminos desea trabajar para n1, solo con enteros se trabaja: "))
        n2 = int(input("Con Cuantos terminos desea trabajar para n2, solo con enteros se trabaja: "))

        x1 = []
        x2 = []

        fx1 = []
        fx2 = []
        a = float(input("Cual es el valor de a: "))
        b = float(input("cual es el valor de b: "))
        h1 = self.TrapezoidWidth(a,b,n1)
        h2 = self.TrapezoidWidth(a,b,n2)
        print("El valor de h1 es: ",h1,"y el valor de h2 es: ",h2)
        x1.append(a)
        x2.append(a)

        for i in range(1,n1):
            x1.append(x1[i-1]+h1)

        for j in range(1,n2):
            x2.append(x2[j-1]+h2)

        x1.append(b)
        x2.append(b)

        f = m.e**(x**2)
        fx1 = [round(s.sympify(f).subs(x,x1[k]),6) for k in range(0,n1+1)]
        fx2 = [round(s.sympify(f).subs(x,x2[l]),6) for l in range(0,n2+1)]

        I1 = self.RichardsonTrape(fx1,n1,x1)
        I2 = self.RichardsonTrape(fx2,n2,x2)
        It = self.balance(I1,I2)
        
        print("Integral aproximada ya con el balanceo",It)
        fireal = lambda x: m.e**(x**2)#np.exp(x**2)
        Irealf = integrate.quad(fireal, a, b)
        print("integral original",Irealf[0])

        E1 = self.Error(Irealf[0],It)
        print("El error es : ",E1," %")
        self.graph(x1,fx1,"Trapecio",x2,fx2)

class Simpson(Help):
    def RichardsonSimpson13(self,fx,n,x,h):
        sfxpares = sfximpares = 0
        for i in range(1,n):
            if i%2==0:
                sfxpares+=fx[i-2]
                print("pares",sfxpares)
            else:
                sfximpares+=fx[i-1]
                print("impares",sfximpares)
        print(fx[0])
        print(fx[n])

        I = (h/3)*(fx[0]+2*sfxpares+4*sfximpares+fx[n])
        return round(I,6)

    def RichardsonSimpson38(self,fx,n,x,h):
        sfx = 0
        for i in range(1,n):
            sfx+=fx[i]
        I = ((3*h)/8)*(fx[0]+3*sfx+fx[n])
        return round(I,6)

    def body13(self):
        x = s.Symbol('x')
        n = int(input("Con Cuantos terminos desea trabajar para n, solo con enteros se trabaja: "))

        X = []

        fx = []
        a = float(input("Cual es el valor de a: "))
        b = float(input("cual es el valor de b: "))
        h = self.TrapezoidWidth(a,b,n)

        print("El valor de h es: ",h)
        X.append(a)

        for i in range(1,n):
            X.append(X[i-1]+h)

        X.append(b)

        f = (1+x**2)**0.5
        fx = [round(s.sympify(f).subs(x,X[k]),6) for k in range(0,n+1)]

        I = self.RichardsonSimpson13(fx,n,X,h)

        print("Integral aproximada: ",I)
        fireal = lambda x: (1+x**2)**0.5#np.exp(x**2)
        Irealf = integrate.quad(fireal, a, b)
        print("integral original",round(Irealf[0],6))

        E1 = self.Error(round(Irealf[0],6),I)
        print("El error es : ",E1," %")
        self.graphSimpson(X,fx,"Simpson 1/3")

    def body38(self):
        x = s.Symbol('x')
        n = int(input("Con Cuantos terminos desea trabajar para n1, solo con enteros se trabaja: "))

        X = []
        fx = []
        a = float(input("Cual es el valor de a: "))
        b = float(input("cual es el valor de b: "))
        h = self.TrapezoidWidth(a,b,n)

        print("El valor de h1 es: ",h)
        X.append(a)

        for i in range(1,n):
            X.append(X[i-1]+h)

        X.append(b)

        f = m.e**(x**2)
        fx = [round(s.sympify(f).subs(x,X[k]),6) for k in range(0,n+1)]

        I = self.RichardsonSimpson38(fx,n,X,h)

        print("Integral aproximada: ",I)
        fireal = lambda x: m.e**(x**2)
        Irealf = integrate.quad(fireal, a, b)
        print("integral original",round(Irealf[0],6))

        E1 = self.Error(round(Irealf[0],6),I)
        print("El error es : ",E1," %")
        self.graphSimpson(X,fx,"Simpson 3/8")

trapezoidal = Trapezoidal()
simpson = Simpson()

numero = 0
print("Extrapolacion de Richardson, e Integracion Numerica")
while (numero != 4):
    numero=int(input("1.Metodo Trapezoidal \n2.Simpson 1/3  \n3.Simpson 3/8  \n4.Salir"))
    if numero==1:
        trapezoidal.body()
    elif numero==2:
        simpson.body13()
    elif numero==3:
        simpson.body38()
    elif numero==4:
        print("Fin ")
