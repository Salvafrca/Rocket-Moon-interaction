program rk
   implicit none
   real*8 t,h, theta, y_1,y_2,y_3,y_4, v, v_, A,mu,Hamiltoniano
   real*8 r,r_,fi, fi_,p_r_,p_fi_ ,pi
   real*8, parameter :: G = 6.67d-11, m = 1.d0 , M_T = 5.9736d24, M_L = 0.07349d24, d_TL = 3.844d8, R_T = 6.37816d6
   real*8, parameter :: omega = 2.6617d-6, R_L = 1.7374d6
   integer n, i, j,corte
   real*8, dimension(1:4, 1:4) :: k
   real*8, external :: f_r_, f_p_r_, f_fi, f_p_fi_

   pi = 4.d0*datan(1.d0)
   t = 0.d0
   h = 0.1d0
   fi = pi/2.d0!Latitud siendo 0 el Polo Norte


   ! Iniciales
   theta = pi/3.39d0
   v = 11200.d0


   !--------------------------------
   !Valores iniciales
   !--------------------------------

   !ORBITA LIGADA NO SE ESCAPA
   ! 1- r = R_T
   ! 1- v = 7903.75d0
   ! 1- theta = 0.d0

   !Casi pero no esscapa
   ! 2- r = R_T
   ! 2- v = 10000.d0
   !2- theta = pi/2.d0

   !SI que escapa
   ! 3- r = R_T
   ! 3-v = 14000.d0
   ! 3-theta = pi/2.d0

   !CHoca con la luna y rebota
   r = R_T
   v = 11200.d0
   theta = pi/3.39d0


   !X,Y, GIF POSICIONES.DAT
   !ENERGIA: GRAFICARLAS EN FUNCION DEL TIEMPO



   ! 0º DEFINICIÓN DE LOS ELEMENTOS DEL VECTOR (REESCALANDOLOS) Y DARLE LOS VALORES INICIALES
   ! LO QUE SE NECESITA PARA TENER CONDICIONES INICIALES ES LO VALORES DE LAS VARIABLES
   ! PARA EL TIEMPO INICIAL, ES DECIR y_1(t_0), y_2(t_0), y_3(t_0), y_4(t_0)
   ! PARA LO QUE SIRVE EL AGORITMO ES PARA QUE DESPUES DE CADA PASO CONOZCA EL VALOR DE
   ! LAS VARIABLES EN EL INSTANTE DE TIEMPO SIGUIENTA DANDO UN PASO DE h
   ! POR TANTO TRAS LA ITERACION CONOCERE y(t_0+h) Y PARA EL SIGUIENTE PASO, CONSIDERARE QUE
   ! y(t_0+h) SON LAS COORDENADAS INICIALES Y POR TANTO CALCULARE y(t_0+h+h) Y ASI SUCESIVAMENTE
   ! x(t) = r(t)*cos(fi(t))
   ! y(t) = y(t)*sen(fi(t))


   ! y_1(t_0)
   ! NECESITO EL VALOR DE LA POSICION INICIAL, A PARTIR DE x(t_0), y(t_0)
   r = R_T !(x*x + y*y)^((1.d0)/(2.d0))
   y_1 = (r)/(d_TL)
   r_ = y_1  ! Radio reescalado

   ! y_2(t_0)
   ! NECESITO EL VALOR DE LA VELOCIDAD INICIAL, A PARTIR DE v_x(t_0), v_y(t_0)
   v_ = v/(d_TL)
   y_2 = v_*cos(theta-fi)
   p_r_ = y_2  ! Momento reescalado

   ! y_3(t_0)
   ! NECESITO EL VALOR DEL ANGULO EN EL INSTANTE INICIAL
   y_3 = fi
   fi_= y_3 ! El angulo no se reescala

   ! y_4(t_0)
   ! NECESITO EL VALOR DEL MOMENTO ANGULAR EN EL INSTANTE INICIAL
   ! PARA ELLO NECESITO EL ANGULO DEL VECTOR VELOCIDAD RESPECTO DE LA HORIZONTAL (theta)
   ! TAMBIÉN NECESITO EL ANGULO DEL VECTOR POSICION RESPECTO DE LA HORIZONTAL (fi)
   p_fi_ = r_*v_*sin(theta-fi)
   y_4 = p_fi_

   ! CON ESTO YA TENDRIA LOS VALORES DE TODAS LAS VARIABLES DEL VECTOR PARA EL INSTANTE INICIAL


   ! DEFINCION DE LAS FUNCIONES y_n PARA PODER SER LLAMADAS DURANTE EL CALCULO NUMERICO

   ! F1 solo usa la varibale p_r_
   !f_r_(p_r_)= p_r_ ! La función para el r reescalado

   ! F2 usa las variables p_fi_, r_, t
   !f_p_r_(p_fi_, r_, t) = (p_fi_*p_fi_)/((r_)^3) - (A)*((1.d0)/(r_*r_) + (mu/((r_)^3))*(r_ - cos(fi- omega*t)))
   ! r_ TILDE!!!!! QUE ES!!!!!

   ! F3 usa las variables p_fi_, r_
   !f_fi(p_fi, r_) = (p_fi_) / ((r_)*(r_))

   ! F4 usa la variable r_, t
   !f_p_fi_(r_, t)= (((-A)*(mu)*(r_)) / ((r_)^3)) * sin(fi - (omega * t))
   open(10,file="Posiciones.dat",status='unknown')
   open(11, file="Energia.dat", status = 'unknown')
   open(12, file="Luna.dat", status = 'unknown')

   !Escribimos posiciones en x e y

   print*, "Iniciales:", r_, fi_,p_fi_,p_r_,v,v_
   A = (G*M_T)/d_TL**3
   mu = M_L/M_T
   corte = 0
   ! TRAS ESTOS PREPARATIVOS PREVIOS YA PUEDO COMENZAR A CALCULAR EL VALOR DE LAS VARIABLES PARA INCREMENTOS DE TIEMPO
   ! CADA ITERACION ME DA LOS VALORES PARA UN VALOR DE TIEMPO DISTINTO.
   ! LA ITERACION 1 ME DA EL VALOR DE y(t+h), LA ITERACION n ME DA EL VALOR DE y(t+n*h)
   do while (t<1000000)

      if (corte == 1000) then
         Hamiltoniano = (p_r_**2)/2.d0 + (p_fi_**2)/(2.d0*r_**2)-A*((1.d0)/r_ + mu/sqrt(r_**2 + 1.d0 -2.d0*r_*cos(fi_-omega*t)))
         write(10,*) r_*cos(fi_), r_*sin(fi_)
         write(11,*) t, Hamiltoniano, Hamiltoniano - omega*p_fi_
         write(12,*) cos(omega*t), sin(omega*t)
         corte = 0
      end if
      corte = corte+1
      ! PARA CALCULAR LAS k(n,m) SE USAN REPETIDAMENTE LAS FUNCIONES f_n Y NOSIEMPRE CON LOS MISMOS VALORES
      ! POR LO QUE PARA SIMPLIFICAR EL CODIGO, SE HAN DEFINIDO PREVIAMENTE COMO SUBFUNCIONES Y POR TANTO
      ! CADA VEZ QUE SEA NECESARIO USARLAS, BASTA CON LLAMARLAS INDICANDO EL VALOR QUE TOMA EN CADA CASO

      ! 2º CALCULO RUNGE-KUTA DE PRIMER ORDEN PARA TODOS LOS ELEMENTOS DEL VECTOR

      k(1,1) = h * f_r_(p_r_)
      k(1,2) = h * f_p_r_(p_fi_, r_,fi_, t,omega,G,M_T,M_L,d_TL)
      k(1,3) = h * f_fi(p_fi_, r_)
      k(1,4) = h * f_p_fi_(r_, fi_,t,omega,G,M_T,M_L,d_TL)


      ! 3º CALCULO DE RUNGE KUTA PARA SEGUNDO ORDEN USANDO LO YA CALCULADO

      k(2,1) = h * f_r_(p_r_ + k(1,2)/2.d0)
      k(2,2) = h * f_p_r_(p_fi_ + k(1,4)/2.d0, r_ + k(1,1)/2.d0,fi_+k(1,3)/2.d0, t + (h)/(2.d0),omega,G,M_T,M_L,d_TL)
      k(2,3) = h * f_fi(p_fi_ + k(1,4)/2.d0, r_ + k(1,1)/2.d0)
      k(2,4) = h * f_p_fi_(r_ + k(1,1)/2.d0,fi_+k(1,3)/2.d0, t + (h)/2.d0,omega,G,M_T,M_L,d_TL)

      ! 4º CALCULO DE RUNGE KUTA PARA TERCER ORDEN USANDO LO YA CALCULADO
      k(3,1) = h * f_r_(p_r_ + k(2,2)/2.d0)
      k(3,2) = h * f_p_r_(p_fi_ + k(2,4)/2.d0, r_ + k(2,1)/2.d0,fi_+k(2,3)/2.d0, t + (h)/(2.d0),omega,G,M_T,M_L,d_TL)
      k(3,3) = h * f_fi(p_fi_ + k(2,4)/2.d0, r_ + k(2,1)/2.d0)
      k(3,4) = h * f_p_fi_(r_ + k(2,1)/2.d0,fi_+k(2,3)/2.d0, t + (h)/2.d0,omega,G,M_T,M_L,d_TL)

      ! 5º CALCULO DE RUNGE KUTA PARA CUARTO ORDEN USANDO LO YA CALCULADO
      k(4,1) = h * f_r_(p_r_ + k(3,2))
      k(4,2) = h * f_p_r_(p_fi_ + k(3,4), r_ + k(3,1),fi_+k(3,3), t + (h),omega,G,M_T,M_L,d_TL)
      k(4,3) = h * f_fi(p_fi_ + k(3,4), r_ + k(3,1))
      k(4,4) = h * f_p_fi_(r_ + k(3,1),fi_ + k(3,3), t + (h),omega,G,M_T,M_L,d_TL)

      ! 6º CALCULO EL VALOR DE CADA VARIABLE PARA EL INSTANTE DE TIEMPO (t + h)

      r_ = r_ + ((1.d0)/(6.d0))*(k(1,1) + (2.d0)*(k(2,1)) + (2.d0)*(k(3,1)) + k(4,1))
      p_r_ = p_r_ + ((1.d0)/(6.d0))*(k(1,2) + (2.d0)*(k(2,2)) + (2.d0)*(k(3,2)) + k(4,2))
      fi_ = fi_ + ((1.d0)/(6.d0))*(k(1,3) + (2.d0)*(k(2,3)) + (2.d0)*(k(3,3)) + k(4,3))
      p_fi_ = p_fi_ + ((1.d0)/(6.d0))*(k(1,4) + (2.d0)*(k(2,4)) + (2.d0)*(k(3,4)) + k(4,4))

      ! 7º AUMENTO EL PASO DEL TIEMPO

      t = t + h

      ! 8º GUARDO EL RESULTADO DE y(t+h) EN CADA ITERACION EN UN FICHERO PARA LUEGO PINTARLO


   end do

   close(10);close(11);close(12)
end program rk

function f_r_(p_r_f) result(resultado)
   implicit none
   real*8 p_r_f,resultado

   resultado = p_r_f

end function

function f_p_r_(p_fi_f, r_f, fi_f,t,omega,G,M_T,M_L,d_TL) result(resultado)
   implicit none
   real*8 p_fi_f,r_f,fi_f, t,resultado, omega
   real*8 r_tilda, A, mu
   real*8 G,M_T,M_L,d_TL
   A = (G*M_T)/d_TL**3
   mu = M_L/M_T
   r_tilda= sqrt(1+ (r_f)**2 - 2.d0*r_f*cos(fi_f -omega*t ))

   resultado = (p_fi_f*p_fi_f)/((r_f)**3) - (A)*((1.d0)/(r_f*r_f) + (mu/((r_tilda)**3))*(r_f - cos(fi_f- omega*t)))

end function

function f_fi(p_fi_f, r_f) result(resultado)
   implicit none
   real*8 p_fi_f,resultado,r_f


   resultado = (p_fi_f) / ((r_f)*(r_f))

end function




function f_p_fi_(r_f,fi_f, t,omega,G,M_T,M_L,d_TL) result(resultado)
   implicit none
   real*8 t,resultado,r_f,omega,fi_f
   real*8 r_tilda, A, mu
   real*8 G,M_T,M_L,d_TL

   r_tilda= sqrt(1+ (r_f)**2 - 2.d0*r_f*cos(fi_f -omega*t ))
   A = (G*M_T)/d_TL**3
   mu = M_L/M_T

   resultado = (((-A)*(mu)*(r_f)) / ((r_tilda)**3)) * sin(fi_f - (omega * t))

end function
