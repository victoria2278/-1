import scipy.integrate #미분방정식을 해석 용도로 쓰이는 함수
import numpy #수치 대입 함수
import matplotlib.pyplot as pyplot #그래프 그리는 함수



def SEIARS_model(y,t,beta,gamma,seta,epsilon,delta1,delta2,omega,zeta): #함수 지정하기

   S,E,I,R,A = y  #변수

   dS_dt = -beta*S*I + epsilon*R #취약군 비율의 변화율
    
   dE_dt = beta*S*I -seta*E #접촉군 비율의 변화율
   
   dI_dt = seta*E - gamma*I - delta1*I #감염군 비율의 변화율

   dR_dt = gamma*I - epsilon*R #회복군 비율의 변화율
   
   dA_dt = omega*I - zeta*R - delta2*I #위중증 환자 비율

   return([dS_dt,dE_dt,dI_dt,dR_dt,dA_dt])
    
S0 = 0.9 #초기 비율
E0 = 0.1
I0 = 0.0
R0 = 0.0 
A0 = 0.0 
beta =0.8 #오미크론 변이의 감염률
gamma = 0.125 #오미크론 변이의 회복률 
seta= 1/4.2 #오미크론 변이의 잠복기의 역수 
epsilon = 1/20 #오미크론 변이의 재감염률
delta1 = 1/1000 #감염군 치사율
delta2 = 1/30 #위중증 환자 치사율
omega = 0.083 #위중증 환자 비율
zeta = 1/70 #중증회복률

t= numpy.linspace(0,100,100000) #0부터 100까지 시간t를 100000으로 쪼갬 
solution = scipy.integrate.odeint(SEIARS_model,[S0,E0,I0,R0,A0],t,args=(beta,gamma,seta,epsilon,delta1,delta2,omega,zeta)) #미분 방정식 해석
solution = numpy.array(solution) #순서를 정함

def SEIRS1_model(y,t,beta1,gamma1,seta1,epsilon1,delta1,delta2,omega,zeta):
    
    S,E,I1,R,A = y

    dS_dt = -beta1*S*I1 + epsilon1*R #미분방정식
    
    dE_dt = beta1*S*I1 -seta1*E

    dI1_dt = seta1*E - gamma1*I1- delta1*I1 

    dR_dt = gamma1*I1 - epsilon1*R 
    
    dA_dt = omega*I1 - zeta*R - delta2*I1

    return([dS_dt,dE_dt,dI1_dt,dR_dt,dA_dt])

S0 = 0.9 
E0 = 0.1
I0 = 0.0
R0 = 0.0 
A0 = 0.0 
beta1 =0.4 #델타 변이의 상수들
gamma1 = 0.1 
seta1= 1/5.8
epsilon1 = 1/100
delta1 = 1/1000
delta2 = 1/30
omega = 0.083
zeta = 1/70

t= numpy.linspace(0,100,100000) #0부터 100까지 시간t를 100000으로 쪼갬 
solution1 = scipy.integrate.odeint(SEIRS1_model,[S0,E0,I0,R0,A0],t,args=(beta1,gamma1,seta1,epsilon1,delta1,delta2,omega,zeta)) 
solution1 = numpy.array(solution1) 

def SEIRS2_model(y,t,beta2,gamma2,seta2,epsilon2,delta1,delta2,omega,zeta):
    
    S,E,I2,R,A = y
# 델타와 오미크론 변이의 중간 지점을 함수로 표현하여 대입함
    dS_dt = -(beta2- (1/(0.2*t+2.5)))*S*I2 + (epsilon2- (1/(0.2*t+25)))*R 
    
    dE_dt = (beta2- (1/(0.2*t+2.5)))*S*I2 -(seta2- (1/(0.2*t+609/40)))*E

    dI2_dt = (seta2- (1/(0.2*t+609/40)))*E - (gamma2- (1/(0.2*t+40)))*I2 - delta1*I2 

    dR_dt = (gamma2- (1/(0.2*t+40)))*I2 - (epsilon2- (1/(0.2*t+25)))*R 
    
    dA_dt = omega*I2 - zeta*R - delta2*I2

    return([dS_dt,dE_dt,dI2_dt,dR_dt,dA_dt])

S0 = 0.9 
E0 = 0.1
I0 = 0.0
R0 = 0.0
A0 = 0.0 
beta2 = 0.8 
gamma2 = 0.125  
seta2 = 10/42 
epsilon2 = 1/20 
delta1 = 1/1000
delta2 = 1/30
omega = 0.083
zeta = 1/70

t= numpy.linspace(0,100,100000) #0부터 100까지 시간t를 100000으로 쪼갬 
solution2 = scipy.integrate.odeint(SEIRS2_model,[S0,E0,I0,R0,A0],t,args=(beta2,gamma2,seta2,epsilon2,delta1,delta2,omega,zeta)) 
solution2 = numpy.array(solution2) 



pyplot.plot(t,solution2[:,0],label="I(t)")
pyplot.plot(t,solution2[:,1],label="I(t)") 
pyplot.plot(t,solution2[:,2],label="I(t)") 
pyplot.plot(t,solution2[:,3],label="I(t)") 
pyplot.plot(t,solution2[:,4],label="I(t)") #오미크론
pyplot.show()
