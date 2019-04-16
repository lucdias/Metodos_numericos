from sympy import *
import matplotlib
import matplotlib.pyplot

#declaracao dos arrays usados para plotar os graficos e para uso dos y's encontrados
plot_x = []
plot_y = []
ys = []
ym = []

#abertura e fechamento dos arquivos
inFile = open('in.txt', 'r')
outFile = open('out.txt','w')

#definicao dos simbolos usados durante todo o programa nas funcoes da biblioteca sympy
y,t = symbols('y t')


#funcao extensa que calcula os devidos y's do metodo de Adams Bashforth
def calcBash(t0,h,expr,order,i):

    if order==2:
        y0 = ys[0+i-order]
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order]
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        return y1 + (h/2)*(3*f1 - f0)
    elif order==3:
        y0 = ys[0+i-order]
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order]
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        y2 = ys[2+i-order]
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        return y2 + (h/12)*(23*f2-16*f1+5*f0)
    elif order==4:
        y0 = ys[0+i-order]
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order]
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        y2 = ys[2+i-order]
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        y3 = ys[3+i-order]
        f3 = expr.subs([(t,t0+(3*h)),(y,y3)])
        return y3 + (h/24)*(55*f3-59*f2+37*f1-9*f0)
    elif order==5:
        #print("-----------")
        y0 = ys[0+i-order]
        #print(y0)
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order]
        #print(y1)
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        y2 = ys[2+i-order]
        #print(y2)
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        y3 = ys[3+i-order]
        #print(y3)
        f3 = expr.subs([(t,t0+(3*h)),(y,y3)])
        y4 = ys[4+i-order]
        #print(y4)
        f4 = expr.subs([(t,t0+(4*h)),(y,y4)])
        #print("-----------")
    
        return y4 + (h)*(((1901/720)*f4)-((1387/360)*f3)+((109/30)*f2)-((637/360)*f1)+((251/720)*f0))
    elif order==6:
        #print("B-----------")
        y0 = ys[0+i-order]
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        #print(y0)
        y1 = ys[1+i-order]
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        #print(y1)
        y2 = ys[2+i-order]
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        #print(y2)
        y3 = ys[3+i-order]
        f3 = expr.subs([(t,t0+(3*h)),(y,y3)])
        #print(y3)
        y4 = ys[4+i-order]
        f4 = expr.subs([(t,t0+(4*h)),(y,y4)])
        #print(y4)
        y5 = ys[5+i-order]
        f5 = expr.subs([(t,t0+(5*h)),(y,y5)])
        #print(y5)
        #print("-----------")
        return y5 + (h*((4277/1440)*f5-(2641/480)*f4+(4991/720)*f3-(3649/720)*f2+(959/480)*f1-(95/288)*f0))
    elif order==7:
        y0 = ys[0+i-order]
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order]
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        y2 = ys[2+i-order]
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        y3 = ys[3+i-order]
        f3 = expr.subs([(t,t0+(3*h)),(y,y3)])
        y4 = ys[4+i-order]
        f4 = expr.subs([(t,t0+(4*h)),(y,y4)])
        y5 = ys[5+i-order]
        f5 = expr.subs([(t,t0+(5*h)),(y,y5)])
        y6 = ys[6+i-order]
        f6 = expr.subs([(t,t0+(6*h)),(y,y6)])
        return y6 + h*((198721/60480)*f6-(18637/2520)*f5+(235183/20160)*f4-(10754/945)*f3+(135713/20160)*f2-(5603/2520)*f1+(19087/6040)*f0)
    elif order==8:
        y0 = ys[0+i-order]
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order]
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        y2 = ys[2+i-order]
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        y3 = ys[3+i-order]
        f3 = expr.subs([(t,t0+(3*h)),(y,y3)])
        y4 = ys[4+i-order]
        f4 = expr.subs([(t+(4*h)),(y,y4)])
        y5 = ys[5+i-order]
        f5 = expr.subs([(t+(5*h)),(y,y5)])
        y6 = ys[6+i-order]
        f6 = expr.subs([(t+(6*h)),(y,y6)])
        y7 = ys[7+i-order]
        f7 = expr.subs([(t,t0+(7*h)),(y,y7)])
        return y7 + h*((16083/4480)*f7-(1152169/120960)*f6+(242653/13440)*f5-(296053/13440)*f4+(2102243/120960)*f3-(115747/13440)*f2+(32863/13440)*f1-(5257/17280)*f0)



#funcao extensa que calcula todos os y's do metodo de Adams Moulton
def calcMult(t0,h,expr,order,i):

    if order==2:
        y0 = ys[0+i-order+1]
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order+1]
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        return y0 + h*((1/2)*f1 + (1/2)*f0)
    elif order==3:
        y0 = ys[0+i-order+1]
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order+1]
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        y2 = ys[2+i-order+1]
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        return y1 + h*((5/12)*f2+(2/3)*f1-(1/12)*f0)
    elif order==4:
        y0 = ys[0+i-order+1]
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order+1]
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        y2 = ys[2+i-order+1]
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        y3 = ys[3+i-order+1]
        f3 = expr.subs([(t,t0+(3*h)),(y,y3)])
        return y2 + ((h/24)*(9*f3+19*f2-5*f1+f0))
    elif order==5:
        #print("-----------")
        y0 = ys[0+i-order+1]
        #print(y0)
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order+1]
        #print(y1)
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        y2 = ys[2+i-order+1]
        #print(y2)
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        y3 = ys[3+i-order+1]
        #print(y3)
        f3 = expr.subs([(t,t0+(3*h)),(y,y3)])
        y4 = ys[4+i-order+1]
        #print(y4)
        f4 = expr.subs([(t,t0+(4*h)),(y,y4)])
        #print("-----------")
    
        return y3 + ((h/720)*(251*f4+646*f3-264*f2+106*f1-19*f0)) 
    elif order==6:
        #print("M-----------")
        y0 = ys[0+i-order+1]
        f0 = expr.subs([(t,t0+(1*h)),(y,y0)])
        #print(y0)
        y1 = ys[1+i-order+1]
        f1 = expr.subs([(t,t0+(2*h)),(y,y1)])
        #print(y1)
        y2 = ys[2+i-order+1]
        f2 = expr.subs([(t,t0+(3*h)),(y,y2)])
        #print(y2)
        y3 = ys[3+i-order+1]
        f3 = expr.subs([(t,t0+(4*h)),(y,y3)])
        #print(y3)
        y4 = ys[4+i-order+1]
        f4 = expr.subs([(t,t0+(5*h)),(y,y4)])
        #print(y4)
        y5 = ys[5+i-order+1]
        f5 = expr.subs([(t,t0+(6*h)),(y,y5)])
        #print(y5)
        #print("-----------")
        return y4 + (h)*((95/288)*f5+(1427/1440)*f4-(113/240)*f3+(241/720)*f2-(173/1440)*f1+(3/160)*f0)
    elif order==7:
        #print("-----------")
        y0 = ys[0+i-order+1]
        #print(y0)
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order+1]
        #print(y1)
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        y2 = ys[2+i-order+1]
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        y3 = ys[3+i-order+1]
        f3 = expr.subs([(t,t0+(3*h)),(y,y3)])
        y4 = ys[4+i-order+1]
        f4 = expr.subs([(t,t0+(4*h)),(y,y4)])
        y5 = ys[5+i-order+1]
        f5 = expr.subs([(t,t0+(5*h)),(y,y5)])
        y6 = ys[6+i-order+1]
        f6 = expr.subs([(t,t0+(6*h)),(y,y6)])
        return y5 + h*((19087/60480)*f6 + (2713/2520)*f5 - (15487/20160)*f4 + (586/945)*f3 - (5737/20160)*f2 + (263/2520)*f1 - (863/60480)*f0)
    elif order==8:
        y0 = ys[0+i-order+1]
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order+1]
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        y2 = ys[2+i-order+1]
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        y3 = ys[3+i-order+1]
        f3 = expr.subs([(t,t0+(3*h)),(y,y3)])
        y4 = ys[4+i-order+1]
        f4 = expr.subs([(t+(4*h)),(y,y4)])
        y5 = ys[5+i-order+1]
        f5 = expr.subs([(t+(5*h)),(y,y5)])
        y6 = ys[6+i-order+1]
        f6 = expr.subs([(t+(6*h)),(y,y6)])
        y7 = ys[7+i-order+1]
        f7 = expr.subs([(t,t0+(7*h)),(y,y7)])
        return y6 + h*((5257/17280)*f7 + (139849/120960)*f6 - (4511/4480)*f5 + (123133/120960)*f4 - (88574/120960)*f3 + (1537/4480)*f2 - (11351/120960)*f1 + (275/24192)*f0)

#funcao extensa que calcula todos os y's do metodo de Formulas inversas da diferenciacao
def calcForm(t0,h,expr,order,i):

    if order==2:
        y0 = ys[0+i-order]
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order]
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        yt = ys[2+i-order]
        ft = expr.subs([(t,t0+(2*h)),(y,yt)])
        return (4/3)*y1 - (1/3)*y0 + (2*h/3)*ft
    elif order==3:
        y0 = ys[0+i-order]
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order]
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        y2 = ys[2+i-order]
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        yt = ys[3+i-order]
        ft = expr.subs([(t,t0+(3*h)),(y,yt)])
        return (18/11)*y2 - (9/11)*y1 + (2/11)*y0 + (6*h/11)*ft
    elif order==4:
        y0 = ys[0+i-order]
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order]
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        y2 = ys[2+i-order]
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        y3 = ys[3+i-order]
        f3 = expr.subs([(t,t0+(3*h)),(y,y3)])
        yt = ys[4+i-order]
        ft = expr.subs([(t,t0+(4*h)),(y,yt)])
        return (48/28)*y3 - (36/25)*y2 + (16/25)*y1 - (3/25)*y0 + (12*h/25)*ft
    elif order==5:
        #print("-----------")
        y0 = ys[0+i-order]
        #print(y0)
        f0 = expr.subs([(t,t0+(0*h)),(y,y0)])
        y1 = ys[1+i-order]
        #print(y1)
        f1 = expr.subs([(t,t0+(1*h)),(y,y1)])
        y2 = ys[2+i-order]
        #print(y2)
        f2 = expr.subs([(t,t0+(2*h)),(y,y2)])
        y3 = ys[3+i-order]
        #print(y3)
        f3 = expr.subs([(t,t0+(3*h)),(y,y3)])
        y4 = ys[4+i-order]
        #print(y4)
        f4 = expr.subs([(t,t0+(4*h)),(y,y4)])
        #print("-----------")
    
        yt = ys[5+i-order]
        ft = expr.subs([(t,t0+(5*h)),(y,yt)])
        return (300/137)*y4 - (300/137)*y3 + (200/137)*y2 - (75/137)*y1 + (12/13)*y0 + (60*h/137)*ft
    elif order==6:
        #print("M-----------")
        y0 = ys[0+i-order]
        f0 = expr.subs([(t,t0+(1*h)),(y,y0)])
        #print(y0)
        y1 = ys[1+i-order]
        f1 = expr.subs([(t,t0+(2*h)),(y,y1)])
        #print(y1)
        y2 = ys[2+i-order]
        f2 = expr.subs([(t,t0+(3*h)),(y,y2)])
        #print(y2)
        y3 = ys[3+i-order]
        f3 = expr.subs([(t,t0+(4*h)),(y,y3)])
        #print(y3)
        y4 = ys[4+i-order]
        f4 = expr.subs([(t,t0+(5*h)),(y,y4)])
        #print(y4)
        y5 = ys[5+i-order]
        f5 = expr.subs([(t,t0+(6*h)),(y,y5)])
        #print(y5)
        #print("-----------")
        yt = ys[6+i-order]
        ft = expr.subs([(t,t0+(6*h)),(y,yt)])
        return (360/147)*y5 - (450/147)*y4 + (400/147)*y3 - (225/147)*y2 + (72/147)*y1 - (10/147)*y0 + (60*h/147)*ft
    
    
#funcao que calcula os pontos do metodo de euler    
def euler(y0,t0,h,n,expr,flag):

    #adicao dos valores iniciais aos devidos arrays de plot e de y's
    ta = t0
    ya = y0
    ys.append(ya)
    plot_x.append(ta)
    plot_y.append(ya)

    #esta flag serve para definir que os metodos que utilizem de euler para a geracao de pontos nao printem nada
    if flag==0:
        outFile.write('Metodo de Euler\n')
        outFile.write('y('+str(t0)+')'+' = '+str(y0)+'\n')
        outFile.write('h = '+str(h)+'\n')
        outFile.write('0 '+str(ya)+'\n')
    for i in range(1,n+1):
        #metodo efetivo, calculado para n passos
        yb = ya + h*expr.subs([(t,ta),(y,ya)])
        ya = yb
        ta = ta + h
        ys.append(yb)
       
        if flag==0:
            plot_x.append(ta)
            plot_y.append(ya)
            outFile.write(str(i)+' '+str(ya)+'\n')


    #devios plots
    if flag==0:
        
        matplotlib.pyplot.xlabel("t")
        matplotlib.pyplot.ylabel("y")     
        matplotlib.pyplot.plot(plot_x, plot_y, 'go')
        matplotlib.pyplot.plot(plot_x, plot_y, 'k:', color='blue')
        matplotlib.pyplot.title('Euler')
        matplotlib.pyplot.show()
        outFile.write("\n")
        outFile.write('\n')
        
    return


#metodo de euler inverso
def eulerInverso(y0,t0,h,n,expr,flag):
    ta = t0
    ya = y0
    ys.append(ya)
    plot_x.append(ta)
    plot_y.append(ya)

    #mesma ideia de euler, por enquanto
    if flag==0:
        outFile.write('Metodo de Euler Inverso\n')
        outFile.write('y('+str(t0)+')'+' = '+str(y0)+'\n')
        outFile.write('h = '+str(h)+'\n')
        outFile.write('0 '+str(ya)+'\n')
    for i in range(1,n+1):
        #a diferenca comeca aqui, estamos prevendo o ponto usato para o calculo da funcao
        yeuler = ya + h*expr.subs([(t,ta),(y,ya)]) #predicted with euler
        yb = ya + h*expr.subs([(t,ta+h),(y,yeuler)])
        ya = yb
        ta = ta + h
        ys.append(yb)
        plot_x.append(ta)
        plot_y.append(ya)
        if flag==0:    
            outFile.write(str(i)+' '+str(ya)+'\n')

    if flag==0:
        matplotlib.pyplot.plot(plot_x, plot_y, 'go')
        matplotlib.pyplot.plot(plot_x, plot_y, 'k:', color='green')
        matplotlib.pyplot.title('Euler Inverso')
        matplotlib.pyplot.show()

        outFile.write('\n')
        outFile.write('\n')
        
    return


#metodo de euler aprimorado
def eulerAprimorado(y0,t0,h,n,expr,flag):
    ta = t0
    ya = y0
    ys.append(ya)
    plot_x.append(ta)
    plot_y.append(ya)

    #mesma ideia do euler inverso e do euler, ate entao
    if flag==0:
        outFile.write('Metodo de Euler Aprimorado\n')
        outFile.write('y('+str(t0)+')'+' = '+str(y0)+'\n')
        outFile.write('h = '+str(h)+'\n')
        outFile.write('0 '+str(ya)+'\n')
    for i in range(1,n+1):
        #calculo da previsao e, logo em seguida, o calculo efetivo do metodo
        yeuler = ya + h*expr.subs([(t,ta),(y,ya)])
        yb = ya + h*expr.subs([(t,ta),(y,ya)])/2 + h*expr.subs([(t,ta+h),(y,yeuler)])/2
        ya = yb
        ta = ta + h
        ys.append(yb)
        plot_x.append(ta)
        plot_y.append(ya)
        if flag==0:   
            outFile.write(str(i)+' '+str(ya)+'\n')

    if flag==0:
        
        matplotlib.pyplot.xlabel("t")
        matplotlib.pyplot.ylabel("y")     
        matplotlib.pyplot.plot(plot_x, plot_y, 'go')
        matplotlib.pyplot.plot(plot_x, plot_y, 'k:', color='orange')
        matplotlib.pyplot.title('Euler Aprimorado')
        matplotlib.pyplot.show()
        outFile.write("\n")
        outFile.write('\n')
        
    return

def rk4(y0,t0,h,n,expr,flag):
    ta = t0
    ya = y0

    #opa, pq flag != -1? pois alem do uso deste metodo cm gerdor de pontos, ele tambem eh usado para prever os pontos explicitos de admas moulton!
    if flag!=-1:
        ys.append(ya)
        plot_x.append(ta)
        plot_y.append(ya)
    if flag==0:
        outFile.write('Metodo de Runge-Kutta de 4ª ordem\n')
        outFile.write('y('+str(t0)+')'+' = '+str(y0)+'\n')
        outFile.write('h = '+str(h)+'\n')
        outFile.write('0 '+str(ya)+'\n')
    for i in range(1,n+1):
        #calculo dos devidos k's dp metodo de runge kutta
        k1 = expr.subs([(t,ta),(y,ya)])
        k2 = expr.subs([(t,ta+(h/2)),(y,ya+(k1*h/2))])
        k3 = expr.subs([(t,ta+(h/2)),(y,ya+(k2*h/2))])
        k4 = expr.subs([(t,ta+(h)),(y,ya+(k3*h))])
        #junto do calculo efetivo do proximo ponto
        yb = ya + (h/6)*(k1+2*k2+2*k3+k4)
        ya = yb
        ta = ta + h
        ys.append(yb)
        if flag!=-1:
            plot_x.append(ta)
            plot_y.append(ya)
        if flag==0:
            outFile.write(str(i)+' '+str(ya)+'\n')
        
    if flag==0:
        
        matplotlib.pyplot.xlabel("t")
        matplotlib.pyplot.ylabel("y")     
        matplotlib.pyplot.plot(plot_x, plot_y, 'go')
        matplotlib.pyplot.plot(plot_x, plot_y, 'k:', color='black')
        matplotlib.pyplot.title('Runge-Kutta de 4ª ordem')
        matplotlib.pyplot.show()
        
        outFile.write('\n')
        outFile.write('\n')
        

    return

#funcao principal para o calculo do metodo de adams bashforth!
def ab(t0,h,n,expr,order,flag):
    ta = t0
    ya = ys[0]
    if flag==2:
        outFile.write('Metodo de Adams-Bashforth de '+str(order)+'ª ordem por Euler\n')
    elif flag==3:
        outFile.write('Metodo de Adams-Bashforth de '+str(order)+'ª ordem por Euler Inverso\n')
    elif flag==4:
        outFile.write('Metodo de Adams-Bashforth de '+str(order)+'ª ordem por Euler Aprimorado\n')
    elif flag==5:
        outFile.write('Metodo de Adams-Bashforth de '+str(order)+'ª ordem por Runge Kutta de 4ª ordem\n')
    elif flag==0:
        outFile.write('Metodo de Adams-Bashforth de '+str(order)+'ª ordem\n')

    outFile.write('y('+str(t0)+') = '+str(ys[0])+'\n')
    outFile.write('h = '+str(h)+'\n')
    outFile.write(str(0)+' '+str(ys[0])+'\n')
    #mostrando os pontos que serao usados
    for i in range(1,order):
        outFile.write(str(i)+' '+str(ys[i])+'\n')
        ta = ta + h
    tk = ta
    ta = t0  

    #tratando para numeros de passos menores que a ordem -> duvida
    if n<order:
        k = n+1
    else:
        k=n+2-order
        
    
    for i in range(1,k):
        #calculo efetivo usando a funcao gigante que esta la em cima
        ya = calcBash(ta,h,expr,order,i+order-1)
        
        ys.append(ya)
        if flag!=-1:
            plot_x.append(ta+tk+h)
            plot_y.append(ya)
            outFile.write(str(i+order-1)+' '+str(ya)+'\n')
        ta = ta + h

    

    #deidos plots    
    matplotlib.pyplot.xlabel("t")
    matplotlib.pyplot.ylabel("y")     
    matplotlib.pyplot.plot(plot_x, plot_y, 'go')
    matplotlib.pyplot.plot(plot_x, plot_y, 'k:', color='black')
    if flag==2:
        matplotlib.pyplot.title('Adams-Bashforth por Euler')
    elif flag==3:
        matplotlib.pyplot.title('Adams-Bashforth por Euler Inverso')
    elif flag==4:
        matplotlib.pyplot.title('Adams-Bashforth por Euler Aprimorado')
    elif flag==5:
        matplotlib.pyplot.title('Adams-Bashforth por Runge-Kutta 4ª ordem')
    elif flag==0:
        matplotlib.pyplot.title('Adams-Bashforth')
    matplotlib.pyplot.show()
    outFile.write('\n')
    outFile.write('\n')
    
    return



def am(t0,h,n,expr,order,flag):
    ta = t0
    ya = ys[0]
    if flag==2:
        outFile.write('Metodo de Adams-Multon de '+str(order)+'ª ordem por Euler\n')
    elif flag==3:
        outFile.write('Metodo de Adams-Multon de '+str(order)+'ª ordem por Euler Inverso\n')
    elif flag==4:
        outFile.write('Metodo de Adams-Multon de '+str(order)+'ª ordem por Euler Aprimorado\n')
    elif flag==5:
        outFile.write('Metodo de Adams-Multon de '+str(order)+'ª ordem por Runge Kutta de 4ª ordem\n')
    elif flag==0:
        outFile.write('Metodo de Adams-Multon de '+str(order)+'ª ordem\n')
    outFile.write('y('+str(t0)+') = '+str(ys[0])+'\n')
    outFile.write('h = '+str(h)+'\n')
    outFile.write(str(0)+' '+str(ys[0])+'\n')

    for i in range(1,order):
        outFile.write(str(i)+' '+str(ys[i])+'\n')
        ta = ta + h
    tk = ta
    ta = t0 

    if n<order:
        k = n+1
    else:
        k=n+2-order
    #mesma ideia do adams bash ate aqui
    for i in range(1,k):
        #outFile.write(ys)
        #outFile.write("rk pegou -> "+str(ys[len(ys)-1]))

        #runge? para prever o ponto explicito e alguns debugs
        rk4(ys[len(ys)-1],ta,h,1,expr,-1)
        #outFile.write("rk gerou -> "+str(ys[len(ys)-1]))
        #outFile.write(ys)


        #calculo efetivo do metodo
        ya = calcMult(ta,h,expr,order,i+order-1)
        ys.pop()
        ys.append(ya)

        plot_x.append(ta+tk+h)
        plot_y.append(ya)
        outFile.write(str(i+order-1)+' '+str(ya)+'\n')
        ta = ta + h

    matplotlib.pyplot.xlabel("t")
    matplotlib.pyplot.ylabel("y")     
    matplotlib.pyplot.plot(plot_x, plot_y, 'go')
    matplotlib.pyplot.plot(plot_x, plot_y, 'k:', color='black')
    if flag==2:
        matplotlib.pyplot.title('Adam Multon por Euler')
    elif flag==3:
        matplotlib.pyplot.title('Adam Multon por Euler Inverso')
    elif flag==4:
        matplotlib.pyplot.title('Adam Multon por Euler Aprimorado')
    elif flag==5:
        matplotlib.pyplot.title('Adam Multon por Runge-Kutta 4ª ordem')
    elif flag==0:
        matplotlib.pyplot.title('Adam Multon')
    matplotlib.pyplot.show()
    outFile.write('\n')
    outFile.write('\n')
    
    return


def fi(t0,h,n,expr,order,flag):
    ta = t0
    ya = ys[0]
    if flag==2:
        outFile.write('Metodo de Formulas inversas da diferenciacao de '+str(order)+'ª ordem por Euler\n')
    elif flag==3:
        outFile.write('Metodo de Formulas inversas da diferenciacao de '+str(order)+'ª ordem por Euler Inverso\n')
    elif flag==4:
        outFile.write('Metodo de Formulas inversas da diferenciacao de '+str(order)+'ª ordem por Euler Aprimorado\n')
    elif flag==5:
        outFile.write('Metodo de Formulas inversas da diferenciacao de '+str(order)+'ª ordem por Runge Kutta de 4ª ordem\n')
    elif flag==0:
        outFile.write('Metodo de Formulas inversas da diferenciacao de '+str(order)+'ª ordem\n')
    outFile.write('y('+str(t0)+') = '+str(ys[0])+'\n')
    outFile.write('h = '+str(h)+'\n')
    outFile.write(str(0)+' '+str(ys[0])+'\n')

    for i in range(1,order-1):
        outFile.write(str(i)+' '+str(ys[i])+'\n')
        ta = ta + h
    tk = ta
    ta = t0 


    if n<order:
        k = n+1
    else:
        k=n+2-order+1

    for i in range(1,k):
        #outFile.write(ys)
        #outFile.write("rk pegou -> "+str(ys[len(ys)-1]))
        rk4(ys[len(ys)-1],t0,h,1,expr,-1)
        #outFile.write("rk gerou -> "+str(ys[len(ys)-1]))
        #outFile.write(ys)

        #quase um adams moulton ate aqui
        ya = calcForm(ta,h,expr,order,i+order-2)
        #pop? Ja que usamos um y previsto para calcular o f explicito, devemos retira-lo do vetor de pontos legitimos
        ys.pop()
        ys.append(ya)

        plot_x.append(ta+tk+h)
        plot_y.append(ya)
        outFile.write(str(i+order-2)+' '+str(ya)+'\n')
        ta = ta + h


    
    matplotlib.pyplot.xlabel("t")
    matplotlib.pyplot.ylabel("y")     
    matplotlib.pyplot.plot(plot_x, plot_y, 'go')
    matplotlib.pyplot.plot(plot_x, plot_y, 'k:', color='black')
    if flag==2:
        matplotlib.pyplot.title('Formulas inversas da diferenciacao por Euler')
    elif flag==3:
        matplotlib.pyplot.title('Formulas inversas da diferenciacao por Euler Inverso')
    elif flag==4:
        matplotlib.pyplot.title('Formulas inversas da diferenciacao por Euler Aprimorado')
    elif flag==5:
        matplotlib.pyplot.title('Formulas inversas da diferenciacao por Runge-Kutta 4ª ordem')
    elif flag==0:
        matplotlib.pyplot.title('Formulas inversas da diferenciacao')
    matplotlib.pyplot.show()
    outFile.write('\n')
    outFile.write('\n')
    
    return
    



def main():

    #leitura das linhas do arquivo
    methods = inFile.readlines()
    
    for i in methods:
        #linha por linha
        m = i.split()
        method = m[0]
        
        
        if m[0] == 'adam_bashforth':
            order = int(m[len(m)-1])
            t0 = float(m[order+1])
            h = float(m[order+2])
            n = int(m[order+3])
            expr = sympify(m[order+4])
            ta = t0
            #gerando a lista inicial de pontos a partir da entrada
            for j in range(1,order+1):
                ys.append(float(m[j]))
                plot_x.append(ta)
                plot_y.append(float(m[j]))
                ta = ta+h
        
        elif m[0] == 'adam_multon':
            order = int(m[len(m)-1])
            t0 = float(m[order+1])
            h = float(m[order+2])
            n = int(m[order+3])
            expr = sympify(m[order+4])
            ta = t0
            #gerando a lista inicial de pontos a partir da entrada
            for j in range(1,order+1):
                ys.append(float(m[j]))
                plot_x.append(ta)
                plot_y.append(float(m[j]))
                ta = ta+h

        elif m[0] == 'formula_inversa':
            order = int(m[len(m)-1])
            t0 = float(m[order])
            h = float(m[order+1])
            n = int(m[order+2])
            expr = sympify(m[order+3])
            ta = t0
            #gerando a lista inicial de pontos a partir da entrada
            for j in range(1,order):
                ys.append(float(m[j]))
                plot_x.append(ta)
                plot_y.append(float(m[j]))
                ta = ta+h
            #print(ys)
            

        
        elif m[0] == 'adam_bashforth_by_euler':
            y0 = float(m[1])
            t0 = float(m[2])
            h = float(m[3])
            n = int(m[4])
            expr = sympify(m[5])
            order = int(m[6])
            #gerando pontos pelo metodo pedido
            euler(y0,t0,h,order-1,expr,1)
            plot_x.clear()
            plot_y.clear()
        
        elif m[0] == 'adam_bashforth_by_euler_inverso':
            y0 = float(m[1])
            t0 = float(m[2])
            h = float(m[3])
            n = int(m[4])
            expr = sympify(m[5])
            order = int(m[6])
            #gerando pontos pelo metodo pedido
            eulerInverso(y0,t0,h,order-1,expr,1)
            plot_x.clear()
            plot_y.clear()

        elif m[0] == 'adam_bashforth_by_euler_aprimorado':
            y0 = float(m[1])
            t0 = float(m[2])
            h = float(m[3])
            n = int(m[4])
            expr = sympify(m[5])
            order = int(m[6])
            #gerando pontos pelo metodo pedido
            eulerAprimorado(y0,t0,h,order-1,expr,1)
            plot_x.clear()
            plot_y.clear()    

        elif m[0] == 'adam_bashforth_by_runge_kutta':
            y0 = float(m[1])
            t0 = float(m[2])
            h = float(m[3])
            n = int(m[4])
            expr = sympify(m[5])
            order = int(m[6])
            #gerando pontos pelo metodo pedido
            rk4(y0,t0,h,order-1,expr,1)
            plot_x.clear()
            plot_y.clear()

        elif m[0] == 'adam_multon_by_euler':
            y0 = float(m[1])
            t0 = float(m[2])
            h = float(m[3])
            n = int(m[4])
            expr = sympify(m[5])
            order = int(m[6])
            #gerando pontos pelo metodo pedido
            euler(y0,t0,h,order-1,expr,1)
            plot_x.clear()
            plot_y.clear()
        
        elif m[0] == 'adam_multon_by_euler_inverso':
            y0 = float(m[1])
            t0 = float(m[2])
            h = float(m[3])
            n = int(m[4])
            expr = sympify(m[5])
            order = int(m[6])
            #gerando pontos pelo metodo pedido
            eulerInverso(y0,t0,h,order-1,expr,1)
            plot_x.clear()
            plot_y.clear()

        elif m[0] == 'adam_multon_by_euler_aprimorado':
            y0 = float(m[1])
            t0 = float(m[2])
            h = float(m[3])
            n = int(m[4])
            expr = sympify(m[5])
            order = int(m[6])
            #gerando pontos pelo metodo pedido
            eulerAprimorado(y0,t0,h,order-1,expr,1)
            plot_x.clear()
            plot_y.clear()    

        elif m[0] == 'adam_multon_by_runge_kutta':
            y0 = float(m[1])
            t0 = float(m[2])
            h = float(m[3])
            n = int(m[4])
            expr = sympify(m[5])
            order = int(m[6])
            #gerando pontos pelo metodo pedido
            rk4(y0,t0,h,order-1,expr,1)
            plot_x.clear()
            plot_y.clear()

        elif m[0] == 'formula_inversa_by_euler':
            y0 = float(m[1])
            t0 = float(m[2])
            h = float(m[3])
            n = int(m[4])
            expr = sympify(m[5])
            order = int(m[6])
            #gerando pontos pelo metodo pedido
            euler(y0,t0,h,order-2,expr,1)
            plot_x.clear()
            plot_y.clear()
        
        elif m[0] == 'formula_inversa_by_euler_inverso':
            y0 = float(m[1])
            t0 = float(m[2])
            h = float(m[3])
            n = int(m[4])
            expr = sympify(m[5])
            order = int(m[6])
            #gerando pontos pelo metodo pedido
            eulerInverso(y0,t0,h,order-2,expr,1)
            plot_x.clear()
            plot_y.clear()

        elif m[0] == 'formula_inversa_by_euler_aprimorado':
            y0 = float(m[1])
            t0 = float(m[2])
            h = float(m[3])
            n = int(m[4])
            expr = sympify(m[5])
            order = int(m[6])
            #gerando pontos pelo metodo pedido
            eulerAprimorado(y0,t0,h,order-2,expr,1)
            plot_x.clear()
            plot_y.clear()    

        elif m[0] == 'formula_inversa_by_runge_kutta':
            y0 = float(m[1])
            t0 = float(m[2])
            h = float(m[3])
            n = int(m[4])
            expr = sympify(m[5])
            order = int(m[6])
            #gerando pontos pelo metodo pedido
            rk4(y0,t0,h,order-2,expr,1)
            plot_x.clear()
            plot_y.clear()

            
        else:
            #tratamento similar para alguns metodos que restaram
            y0 = float(m[1])
            t0 = float(m[2])
            h = float(m[3])
            n = int(m[4])
            expr = sympify(m[5])
            order = -1
           
        #DEBUGGGGGG
        #print("\nEsse eh o metodo = "+method)
        #if method!='adam_bashforth' and method!='adam_multon' and method!='formula_inversa':
        #    print("Esse eh o y0 = "+str(y0))
        #print("Esse eh o t0 = "+str(t0))
        #print("Esse eh o passo = "+str(h))
        #print("Esse eh o numero de passos = "+str(n))
        #print("Essa eh a expressao = "+str(expr))
        #print("Essa eh a ordem = "+str(order))

        #chamada efetiva para eles
        if method=='euler':
            euler(y0,t0,h,n,expr,0)
            ys.clear()
        elif method=='euler_inverso':
            eulerInverso(y0,t0,h,n,expr,0)
            ys.clear()
        elif method=='euler_aprimorado':
            eulerAprimorado(y0,t0,h,n,expr,0)
            ys.clear()
        elif method=='runge_kutta':
            rk4(y0,t0,h,n,expr,0)
            ys.clear()
        elif method=='adam_bashforth':
            ab(t0,h,n,expr,order,0)
            ys.clear()
        elif method=='adam_bashforth_by_euler':
            ab(t0,h,n,expr,order,2)
            ys.clear()
        elif method=='adam_bashforth_by_euler_inverso':
            ab(t0,h,n,expr,order,3)
            ys.clear()
        elif method=='adam_bashforth_by_euler_aprimorado':
            ab(t0,h,n,expr,order,4)
            ys.clear()
        elif method=='adam_bashforth_by_runge_kutta':
            ab(t0,h,n,expr,order,5)
            ys.clear()
        elif method=='adam_multon':
            am(t0,h,n,expr,order,0)
            ys.clear()
        elif method=='adam_multon_by_euler':
            am(t0,h,n,expr,order,2)
            ys.clear()
        elif method=='adam_multon_by_euler_inverso':
            am(t0,h,n,expr,order,3)
            ys.clear()
        elif method=='adam_multon_by_euler_aprimorado':
            am(t0,h,n,expr,order,4)
            ys.clear()
        elif method=='adam_multon_by_runge_kutta':
            am(t0,h,n,expr,order,5)
            ys.clear()

        elif method=='formula_inversa':
            fi(t0,h,n,expr,order,0)
            ys.clear()
        elif method=='formula_inversa_by_euler':
            fi(t0,h,n,expr,order,2)
            ys.clear()
        elif method=='formula_inversa_by_euler_inverso':
            fi(t0,h,n,expr,order,3)
            ys.clear()
        elif method=='formula_inversa_by_euler_aprimorado':
            fi(t0,h,n,expr,order,4)
            ys.clear()
        elif method=='formula_inversa_by_runge_kutta':
            fi(t0,h,n,expr,order,5)
            ys.clear()

        plot_x.clear()
        plot_y.clear()
        
        
        

    inFile.close()



main()    
