#! /usr/bin/env python
import rospy
from std_msgs.msg import String, Float32, Int32,Int32MultiArray,MultiArrayLayout,MultiArrayDimension, Float32MultiArray
from nav_msgs.msg import Odometry
from tf.transformations import euler_from_quaternion
from geometry_msgs.msg import Point, Twist
from sensor_msgs.msg import LaserScan
from collections import namedtuple
import math
import numpy as np


def pol2cart(rho, phi):
    N = np.size(rho,0)
    x = rho[0:N] * np.cos(phi[0:N])
    y = rho[0:N] * np.sin(phi[0:N])

    # Stack row on row
    XY =np.vstack((x,y))

    return XY

def fit_line(XY):
    # Number of columns
    length = np.size(XY,1)

    # Centroids
    xc = sum(XY[0])/length
    yc = sum(XY[1])/length

    # Increments
    dX = (XY[0] - xc)
    dY = (XY[1] - yc)

    # Alpha
    num = -2*sum(np.multiply(dX,dY))
    denom = sum(np.multiply(dY, dY) - np.multiply(dX,dX))
    alpha = math.atan2(num, denom)/2

    # r
    r = xc*math.cos(alpha)+yc*math.sin(alpha)

    # Eliminate negative radii
    
    if (r<0):
        alpha = alpha + math.pi
        if (alpha>math.pi):
            alpha = alpha - 2*math.pi
            r = -r

    # Return r and alpha parameters of a line
    return r, alpha

def find_split_point(XY, r, alpha, parameters): # XY is temporary set, NOT GLOBAL !!!!
    # Calculating distances of all points from a line

    N = np.size(XY,1)

    # Provjera koliko segment ima tacaka, ako ih je manje od 10, ne zelimo da splitujemo taj segment
    if (N<=10):
        split = 0
        SplitPointIdx = 0
        return SplitPointIdx, split, N

    cosA = math.cos(alpha)
    sinA = math.sin(alpha)
    d = np.zeros((1,N)) #1xN

    xcosA = XY[0]*cosA 
    ysinA = XY[1]*sinA
    d = xcosA + ysinA - r
    # I need absolute values
    dabs = map(abs, d)
    
    # Find max
    D = np.max(dabs)
    # Check if it is far enough
    if (D>parameters.D):
        # Find index of that point
        split = 1
        result = np.where(dabs==D)
        
        SplitPointIdx = result[0][0]

    else:
        split = 0
        SplitPointIdx = 0



    #print SplitPointIdx, split
    return SplitPointIdx, split, N

def split_lines_r(XY, startIdx, endIdx, parameters): # This returns rs, alphas and start and end points of segments
    global r
    global alpha
    global Idx
    global isk
    global Outliers
    global outl_counter
    N = endIdx - startIdx + 1# Number of elements that set contains

    rac, alpaca = fit_line(XY[:,startIdx:endIdx+1]) # Fit a line to that set

    SplitPointIdx, split, N_temp= find_split_point(XY[:,startIdx:endIdx+1], rac, alpaca, parameters) # Find a point for splitting

    if (split):
        #Prvo provjeriti da li je SplitPointIdx pocetni ili krajnji indeks, tj da li je 0 ili N-1, ukoliko jeste, ne splitujemo tu, nego pozivamo split_lines_r
	if (SplitPointIdx==0):
	    Outliers.append(SplitPointIdx+startIdx) # Dodala sam taj indeks u listu outliera
            #print 'SLR, SP = 0, Autlajeri su:', Outliers
	    outl_counter = outl_counter + 1
	    #Potrebno je pozvati sada ovu funkciju ali bez ovog indeksa koji mi je vracen
	    split_lines_r(XY, startIdx+1, endIdx,parameters)
	elif (SplitPointIdx==(N_temp-1)):
	    Outliers.append(SplitPointIdx+startIdx)
            #print 'SLR, SP = 1, Autlajeri su:', Outliers
	    outl_counter = outl_counter + 1
            #Potrebno je pozvati sada ovu funkciju ali bez ovog indeksa koji mi je vracen
	    split_lines_r(XY, startIdx, endIdx-1,parameters)
	else:	
            split_lines_r(XY, startIdx, SplitPointIdx+startIdx, parameters)
            split_lines_r(XY, SplitPointIdx+startIdx, endIdx, parameters)

    else: # Line cannot be split, save start and end point
        r.insert(isk,rac)
        alpha.insert(isk, alpaca)
        Idx.insert(isk,[startIdx, endIdx])
        isk  =  isk + 1


    return

def merge(XY, rho, phi, r, alpha, Idx, coord_Idx):
    s = [r[0], alpha[0]]
    startIdx = Idx[0][0]
    lastEndIdx = Idx[0][1] # Ovo je indeks posljednje tacke prvog segmenta koji sam uzela
    pom1 = coord_Idx[0] # Ovdje su x,y koordinate prve tacke prvog segmenta
    pom2 = coord_Idx[1] # Ovdje su x,y koordinate posljednje tacke prvog segmenta 
    N = len(r)
	
    rOut = []
    alphaOut = []
    IdxOut = []
    coord_IdxOut = []
    for i in range(1,N): 	   
        endIdx = Idx[i][1] # Indeks posljednje tacke drugog segmenta
        x0 = pom1[0] 
	y0 = pom1[1]
        x1 = pom2[0]
	y1 = pom2[1]
        x2 = coord_Idx[i*2][0]
	y2 = coord_Idx[i*2][1]
	x3 = coord_Idx[(i*2)+1][0]
	y3 = coord_Idx[(i*2)+1][1]
		
	area1 = abs(x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1))
	area2 = abs(x0*(y2-y3)+x1*(y3-y0)+x3*(y0-y2))
	
	# Sada ide provjera
	if ((area1<=0.2) and (area2<=0.3)):
        #if ((area1+area2)<1):
	    # U ovom slucaju smatram da su dva segmenta koje posmatram dovoljno
	    # kolinearni i da je potrebno da fitujem liniju izmedju 
	    # pocetne tacke prvog segmenta i posljednje tacke drugog segmenta
					
	    r_t, alpha_t = fit_line(XY[:, startIdx:endIdx+1])
		
            s = [r_t, alpha_t]
	    # Sada je potrebno samo da nadjem koordinate krajnjih tacaka ovog segmenta
	    indeks = [startIdx,endIdx]
	    x_m = r_t*math.cos(alpha_t)
	    y_m = r_t*math.sin(alpha_t)
	    X = []
	    Y = []
	    for j in range (0,2):
	        a = indeks[j]
		x_p = (abs(r_t)+(rho[a]*math.cos(phi[a]-alpha_t)-abs(r_t)))*math.cos(alpha_t)
	        y_p = (abs(r_t)+(rho[a]*math.cos(phi[a]-alpha_t)-abs(r_t)))*math.sin(alpha_t)
	        x = XY[0][a]
                y = XY[1][a]			
	        dis_x = abs(x_p - x)
	        dis_y = abs(y_p - y)
		if (x_p<=x):
                    x_end = x_m + dis_x
                    X.append(x_end)
                    if (y_p>=y):
                        y_end = y_m - dis_y
                        Y.append(y_end)
                    elif (y_p<=y):					
                        y_end = y_m + dis_y
                        Y.append(y_end)
                    		
	        elif (x_p>=x):
		    x_end = x_m - dis_x
                    X.append(x_end)
                    if (y_p>=y):
                        y_end = y_m - dis_y
                        Y.append(y_end)
                    elif (y_p<=y):					
                        y_end = y_m + dis_y
                        Y.append(y_end)
                    
	        elif (x_m==x and y_m==y):
		    x_end = x
                    X.append(x_end)
		    y_end = y
                    Y.append(y_end)
	    pom1 = [X[0],Y[0]]
	    pom2 = [X[1],Y[1]]
			
        else: # Nije potrebno mergeovanje i samo treba da upisem informacije o segmentu
            rOut.append(s[0])
	    alphaOut.append(s[1])
	    IdxOut.append([startIdx, lastEndIdx])
	    coord_IdxOut.append(pom1)
	    coord_IdxOut.append(pom2)
			
	    # Potrebno se pomjeriti na sljedeci segment
	    s = [r[i],alpha[i]]
	    startIdx = Idx[i][0]
	    pom1 = coord_Idx[i*2] # Ovdje su x,y koordinate prve tacke prvog segmenta
            pom2 = coord_Idx[(i*2)+1] # Ovdje su x,y koordinate posljednje tacke prvog segmenta
	lastEndIdx = endIdx

    rOut.append(s[0])
    alphaOut.append(s[1])
    IdxOut.append([startIdx, lastEndIdx])
    coord_IdxOut.append(pom1)
    coord_IdxOut.append(pom2)		
	# Mozda ce biti potrebno da se doda posljednji segment, a mozda i n
    return rOut, alphaOut, IdxOut, coord_IdxOut

def SAM(rho_temp, phi_temp):

     print('=======================BEGINNING=======================')

     # Parameter strucure:
     parameters = namedtuple("parameters", "D")
     # D is threshold value for splitting
     params = parameters(D = 0.09) # This is parameter for splitting

     # These are final lists, without inf, -inf or nan
     rho = []
     phi = []
     j = 0 # Counter for final lists
     n = 0 # Counter for number of deleted points (because they are inf, -inf, or nan)
     for i in range(0, len(rho_temp)):
    	if ((rho_temp[i] != float('inf')) and (rho_temp[i] != float('-inf')) and (rho_temp[i] != float('nan'))):
    		rho.append(rho_temp[i])
    		phi.append(phi_temp[i])
    		j +=1
    	else:
    		n +=1

     print 'Length of rho:', len(rho)
     print 'Min of rho:', min(rho)
     print 'Max of rho:', max(rho)
     print 'Length of phi:', len(phi)
     print 'Number of deleted points:', n
     print 'Threshold value for splitting:', params.D
     XY  = pol2cart(rho, phi) # Starting set of points

     print 'Number of columns of our starting set:', np.size(XY,1)

     print('=======================SPLITTING=======================')

     global r
     global alpha
     global Idx
     global isk
     global Outliers
     global outl_counter
     r = []
     alpha = []
     Idx = []
     isk = 0
     Outliers = []
     outl_counter = 0

     split_lines_r(XY, 0, np.size(XY,1)-1, params)

     print 'Set of rs: ', r
     print 'Set of alphas: ', alpha
     print 'Indices are: ', Idx
     print 'Number of segments: ', isk

     #plt.plot(XY[[0],:],XY[[1],:], 'ro')
     #print 'Outlajer: ', Outliers
     #print 'Broj outlajera: ', outl_counter
     
     Outliers.sort()
     clones  = [] # Lista indeksa u listi Outliersa onih brojeva koji se po drugi put pojavljuju u listi Outliers
     for i in range(0, len(Outliers)-1):
         # Provjeri da li je trenutni element jednak sljedecem elementu, ako nije, ubaci ga u clones
	 if (Outliers[i]==Outliers[i+1]):
	     clones.append(i+1)
     for i in range(0, len(clones)):
         a = clones[i]
         Outliers.pop(a)
     print 'Autlajeri sad: ', Outliers
     print 'Broj tacaka u Outliers: ', len(Outliers)

     # Sad je potrebno da napravim podsetove gdje cu za svaki da pozovem split_lines_r
     # Najbolje bi bilo da napravim listu listi, i onda da idem kroz listu Outliersa, i onda da imam 2 brojaca, start i end, jedan ce da mi pamti od kog indeksa u listi Outliers
     # sam krenula , a drugi ce da prati tok kretanja
     sets = [] # Ovo je ta prazna lista listi 
     start = 0
     end = 0 
     br = 0
     for i in range(0, len(Outliers)-1):
         current = Outliers[i]
	 next = Outliers[i+1]
	 if (next == (current+1)):
             if (i == (len(Outliers)-2)):
	         end = i+1
                 sets.insert(br, Outliers[start:end+1])
	         br = br+1
             else:
	         end = i+1
	 elif (next!=(current+1)):
             sets.insert(br, Outliers[start:end+1])
	     br = br+1
	     start = end + 1
             end = end + 1
     # Ovo valjda radi, ako ne, POPRAVITI, ali svakako provjeriti svaki ovaj korak
     # Provjer:
     print 'Da li mi ovo izbacuje broj ukupnih setova ili samo broj elemenata jednog seta', len(sets)
     print 'A lista setova: ', sets
     Outliers = []
     outl_counter = 0
     for i in range(0, len(sets)):
         set = sets[i]
         start = set[0]
         end = set[-1]
         if (len(set)>1):
	     split_lines_r(XY, start, end, params)

     Outliers.sort()
     clones  = [] # Lista indeksa u listi Outliersa onih brojeva koji se po drugi put pojavljuju u listi Outliers
     for i in range(0, len(Outliers)-1):
         # Provjeri da li je trenutni element jednak sljedecem elementu, ako nije, ubaci ga u clones
	 if (Outliers[i]==Outliers[i+1]):
	     clones.append(i+1)
     for i in range(0, len(clones)):
         a = clones[i]
         Outliers.pop(a)
     sets = [] # Ovo je ta prazna lista listi 
     start = 0
     end = 0 
     br = 0
     for i in range(0, len(Outliers)-1):
         current = Outliers[i]
	 next = Outliers[i+1]
	 if (next == (current+1)):
             if (i == (len(Outliers)-2)):
	         end = i+1
                 sets.insert(br, Outliers[start:end+1])
	         br = br+1
             else:
	         end = i+1
	 elif (next!=(current+1)):
             sets.insert(br, Outliers[start:end+1])
	     br = br+1
	     start = end + 1
             end = end + 1
     # Ovo valjda radi, ako ne, POPRAVITI, ali svakako provjeriti svaki ovaj korak
     # Provjer:
     #print 'Da li mi ovo izbacuje broj ukupnih setova ili samo broj elemenata jednog seta', len(sets)
     #print 'A lista setova: ', sets
     Outliers = []
     outl_counter = 0
     for i in range(0, len(sets)):
         set = sets[i]
         start = set[0]
         end = set[-1]
         if (len(set)>1):
	     split_lines_r(XY, start, end, params)
     #print 'Autlajeri sad: ', Outliers
     #print 'Broj tacaka u Outliers: ', len(Outliers)
     print 'r je sada: ', r
     print 'Broj linija je sada: ', len(r)
     '''
     set = sets[0]
     start = set[0]
     print start
     end = set[-1]
     print end
     print 'XY', XY[:, start:end+1]'''
     
     '''print 'Autlajeri sad: ', Outliers
     print 'Broj autlajera je: ', len(Outliers)
     print 'Alpha prije sort ', alpha
     print 'Rovi prije sortianja', r
     print 'Indices of segments: ', Idx'''
     Idx_clone = Idx
     Idx, r = (list(t) for t in zip(*sorted(zip(Idx,r))))
     Idx_clone, alpha = (list(t) for t in zip(*sorted(zip(Idx_clone,alpha))))
    

     print 'Alpha nakon sort: ', alpha
     print 'R nakon sortiranja ', r
     print 'Nakon sortiranja', Idx
     # Ovo do sada je sve dobro
     coord_Idx = []
     
     for i in range(0, len(r)):
        X = []
	Y = []
	# This point lies on the center of a segment
	idx = Idx[i] # Indices of first and last point of ith segment
        # Koordinate tacke koja se nalazi na liniji sa parametrima ovog segmenta pod normalnim uglom od robota
	x_m = abs(r[i])*math.cos(alpha[i])
	y_m = abs(r[i])*math.sin(alpha[i])
	for j in range(0,2):
            a = idx[j]

	    # Moram se vratiti u XY set
            # 
	    x_p = (abs(r[i])+(rho[a]*math.cos(phi[a]-alpha[i])-abs(r[i])))*math.cos(alpha[i])
	    y_p = (abs(r[i])+(rho[a]*math.cos(phi[a]-alpha[i])-abs(r[i])))*math.sin(alpha[i])
	    x = XY[0][a]
            y = XY[1][a]
			
	    dis_x = abs(x_p - x)
	    dis_y = abs(y_p - y)
   
	    if (x_p<=x):
                x_end = x_m + dis_x

                X.append(x_end)
                
                if (y_p>=y):
                    y_end = y_m	- dis_y

                    Y.append(y_end)
                    
                elif (y_p<=y):					
                    y_end = y_m	+ dis_y

                    Y.append(y_end)
                    		
	    elif (x_p>=x):
		x_end = x_m - dis_x

                X.append(x_end)
                
                if (y_p>=y):

                    y_end = y_m	- dis_y
                    Y.append(y_end)

                elif (y_p<=y):					
                    y_end = y_m	+ dis_y

                    Y.append(y_end)

	    elif (x_m==x and y_m==y):
		x_end = x
                
		X.append(x_end)
		y_end = y
                
                Y.append(y_end)
        coord_Idx.append([X[0],Y[0]])
        coord_Idx.append([X[1],Y[1]]) 
	#if (X!=[]):
        #   plt.plot([X[0],X[1]],[Y[0],Y[1]], 'bo-')
     #plt.show()
     print('=======================MERGING=======================')

     rOut, alphaOut, IdxOut, coord_IdxOut = merge(XY, rho, phi, r, alpha, Idx, coord_Idx)
     print 'rOut: ', rOut
     print 'alphaOut: ', alphaOut
     print 'IdxOut: ', IdxOut
     print 'Number of segments: ', np.size(rOut,0)
     r = rOut
     alpha = alphaOut
     Idx = IdxOut
     
     return alpha, r
# Bice potrebno da koristimo rqt_plot - istraziti na ROSWiki kako se to koristimo
# Sto se tice rqt_plot kada napisemo nasu glavnu funkciju koja se poziva na svake 0.2 sekunde, samo je potrebno u njoj kreirati Publisher,
# I taj Publisher ce na kraju da salje krajnje podatke (ZNACI PODATKE O POZICIJI ROBOTA KOJU SMO ODREDILI PREKO EKF i ONE PODATKE O POZICIJI KOJE SMO DOBILI SENZOROM 
# kako bismo mogli da ih poredimo) na neki topic, i samo cemo pozivanjem rqt_plot /paimetogtopica dobiti iscrtano to sto Nikola valjda trazi
def measurement_prediction(x, i): # i sluzi da bismo znali koju kolonu da uzimamo iz world_map [dole ima neka for petlja]
    # PROVERITI DA LI JE TO i NA DOBROJ POZICIJI
    # Ova funkcija sluzi za prebacivanje globalnih koordinata u lokalni koordinatni sistem robota, pretpostavljajuci da je on na poziciji 
    # koju smo procenili pomocu odometrije
    global world_map
    h1 = world_map[0][i]-x[2]
    h2 = world_map[1][i] - (x[0]*math.cos(world_map[0][i]) + x[1]*math.sin(world_map[0][i]))
    z_p = [h1,h2] 
    # Jakobijan ove dvije funkcije gore
    H_p= [[0,0,-1],[-math.cos(world_map[0][i]),-math.sin(world_map[0][i]),0]]
	
    # Ovdje mozda provjeriti znakove zbog ovih uglova i r-ova
    return z_p, H_p
	
def state_prediction(x, u, b, Q, P):
    # Ovdje izracunamo predikciju narednog stanja i predikciju kovarijacione matrice
    x_priori = x + np.array([((u[0]+u[1])/2) * math.cos(x[2] + ((u[1]-u[0])/(2*b))), ((u[0]+u[1])/2) * math.sin(x[2] + ((u[1]-u[0])/(2*b))), (u[1]-u[0])/b])
	
    # Jakobijan Fx
    df1_dx = 1
    df1_dy = 0
    df1_dteta = -math.sin(x[2]- (u[0] - u[1])/(2*b))*(u[0]/2 + u[1]/2)
    df2_dx = 0
    df2_dy = 1
    df2_dteta = math.cos(x[2] - (u[0] - u[1])/(2*b))*(u[0]/2 + u[1]/2)
    df3_dx = 0
    df3_dy = 0
    df3_dteta = 1
    F_x = np.array([[df1_dx, df1_dy, df1_dteta],[df2_dx, df2_dy, df2_dteta],[df3_dx, df3_dy, df3_dteta]])
	
    # Jakobijan Fu
    df1_du1 = math.cos(x[2] - (u[0] - u[1])/(2*b))/2 + (math.sin(x[2] - (u[0] - u[1])/(2*b))*(u[0]/2 + u[1]/2))/(2*b)
    df1_du2 = math.cos(x[2] - (u[0] - u[1])/(2*b))/2 - (math.sin(x[2] - (u[0] - u[1])/(2*b))*(u[0]/2 + u[1]/2))/(2*b)
    df2_du1 = math.sin(x[2] - (u[0] - u[1])/(2*b))/2 - (math.cos(x[2] - (u[0] - u[1])/(2*b))*(u[0]/2 + u[1]/2))/(2*b)
    df2_du2 = math.sin(x[2] - (u[0] - u[1])/(2*b))/2 + (math.cos(x[2] - (u[0] - u[1])/(2*b))*(u[0]/2 + u[1]/2))/(2*b)
    df3_du1 = -1/b
    df3_du2 = 1/b
    F_u = np.array([[df1_du1,df1_du2],[df2_du1,df2_du2],[df3_du1,df3_du2]])
	
    # Predikcija kovarijacione matrice
    #P_priori = F_x @ P @ F_x.T + F_u @ Q @ F_u.T
    P_priori = np.linalg.multi_dot([F_x,P,F_x.T]) + np.linalg.multi_dot([F_u,Q,F_u.T])
    return x_priori, P_priori

	
def glavna(rho_temp, phi_temp):
    print 'GLAVNA'
    
    global v
    print 'Brzina v: ', v
    global w
    print 'Brzina w: ', w
    global vl
    print 'Brzina vl: ', vl
    global vr
    print 'Brzina vr: ', vr
    global x
    global world_map
    global first
    global P

    if (first==1):
        alpha, r = SAM(rho_temp, phi_temp)
        world_map = np.vstack((alpha,r))
	first = 0
    print 'Referentna mapa: ', world_map
    global R
	
    sigma = 0.05/3
    R = np.array([[sigma**2,0],[0,sigma**2]])
    print 'Matrica R: ', R
    # Potrebno je da znam k za matricu Q i l, gdje je l rastojanje izmedju tockova
    # Ovu funkciju treba da pozivamo na svake 0.2 sekunde samo ne znam kako
    k = 1 # OVDE JE LUPLJENA VRIJEDNOST BY SRKI - OVO CE TREBATI DA SE SMANJI
    b = 0.16 # rastojanje od jednog do drugog tocka - 160mm GORE NA NEKOM MESTU PISE DA JE L POLOVINA RASTOJANJA, A OVDE DA JE RASTOJANJE??? DA NIJE GRESKA?!
    # Mislim da nije greska posto se trazi rastojanje u ovom dijelu, a kod onog proracuna se trazi polu rastojanje
    # Ovdje cu da kreiram publisher, koji ce na neki topic CRTANJE, da salje podatke, koje cemo iscrtavati preko rqt_plot naredbe
    pub = rospy.Publisher("/crtanje", Float32MultiArray, queue_size = 1) # Treba provjeriti tip podataka da li je ok
	
    # Sada cu da zamislim da imam sve podatke koji su potrebni, znaci ona MAPA i podaci sa odoma, i implementiracu Kalmanov filtar kao da mi je to sve poznato
    delta_t = 0.2 # Valjda treba ovoliko da iznosi
    delta_s_l = vl*delta_t
    print 'DeltaSL: ', delta_s_l
    delta_s_r = vr*delta_t
    print 'DeltaSR: ', delta_s_r
	
    u = [delta_s_l, delta_s_r]
	
    # Prvo ucitam matricu Q
    Q = k*np.diag(map(abs,u)) # Valjda moze ovako
    print 'Matrica Q: ', Q
    # Ovaj ovdje korak se zasniva na prethodnoj procjeni koju je izbacio Kalmanov filtar i odometriji - 1. korak: ACT
    print '========1.korak - ACT========'
    x_priori, P_priori = state_prediction(x, u, b, Q, P)
    print 'x_priori: ', x_priori
    print 'P_priori:', P_priori
    # 2. korak: Opservacija (SEE)
    print '========2.korak - SEE========'
    # Ovdje ocitavamo mjerenja senzora u sadasnjem trenutku, a to podrazumijeva transformaciju globalnih obiljezja mape u lokalni koordinatni sistem 
    # Ovdje je potrebno kopirati onaj kod iz prethodnog zadatka da bi se dobili podaci o mapi
    # !!!
    alpha, r = SAM(rho_temp, phi_temp)
    Z = np.vstack((alpha,r)) # Prva vrsta alpha, druga vrsta r OVDE JE RANIJE BILO M UMESTO Z
    print 'Mjerenja: ', Z
	
	
    # Specifikacije za Lidar senzor:
    # http://emanual.robotis.com/docs/en/platform/turtlebot3/appendix_lds_01/

	
    # Ovdje dalje SMATRAM DA IMAMO TU MATRICU R... 
	
	
    # Mozda je bolje 3 i 4 korak spojiti, da ne bi imali dupli proracun, ali cu zasad da ovo napisem odvojeno
	
    # 3. korak: Predikcija mjerenja 
    print '========3.korak - Predikcija mjerenja========'
    nMapEntries = np.size(world_map,1) # Ovdje trazim broj kolona, tj. broj segmenata koji nam kod iz treceg zadatka daje
    print 'Broj linija u referentnoj mapi: ',nMapEntries
    z_priori = np.zeros((2,nMapEntries)) # Neka bude vektor vrsta
    H = np.zeros((2,3,nMapEntries)) # 3D matrica, jer svaki ovaj z_priori ima svoju H matricu
    for i in range(0, nMapEntries):
	z_priori[:,i], H[:,:,i] = measurement_prediction(x_priori, i)
	# z_priori su parametri linija MAPE predstavljeni u lokalnom koordinatnom sistemu robota na trenutno estimiranoj poziciji!
    print 'Z_priori: ', z_priori   
    # 4. korak: Uparivanje - za uparivanje koristimo razliku izmedju svih parametara linija dobijenih pomocu mjerenja i dobijenih iz
	
    # globalne mape konvertovanjem u lokalni koordinatni sistem, tzv. inovacije
    # Da bismo mogli upariti obiljezja, Mahalonobisova distanca mora biti manja od neke validacione granice
   

    print '========4.korak - Uparivanje mjerenja========'
    nMeasurements = np.size(Z,1)
    d = np.zeros((nMeasurements, nMapEntries))
    v_ino = np.zeros((2, nMeasurements * nMapEntries)) 			#2,n*n dimenzija # <-- JEL OVO TREBALO DA BUDE ZAKOMENTARISANO?
    for i in range(0, nMeasurements): # Iteriramo kroz mjerenja 
        for j in range(0, nMapEntries): # Iteriramo kroz sva prediktovanja mjerenja
	    v_ino[:, j+i*nMapEntries] = Z[:,i] - z_priori[:,j] # inovacija
	    #W = H[:,:,j] @ P_priori @ H[:,:,j].T + R
            W = np.linalg.multi_dot([H[:,:,j],P_priori,H[:,:,j].T]) + R 
	    #d[i][j] = v[:, j+i*nMapEntries].T @ np.linalg.inv(W) @ v[:, j+i*nMapEntries]
            d[i][j] = np.linalg.multi_dot([v_ino[:, j+i*nMapEntries].T,np.linalg.inv(W),v_ino[:, j+i*nMapEntries]])
	    # Kolone matrice d idu po z_priori a vrste po mjerenjima Z
    print 'v_ino prije racunanja distance', v_ino
    print 'Matrica d: ', d
    # Kada smo izracunalu ove inovacije i Mahalonobisove distance, potrebno je u svakoj vrsti pronaci
    # Minimum, i ako je taj minimum manje od zadatog g, onda ta mjerenja uparimo, sacuvamo valjda to v...
    # I nesto se dalje radi, uglavnom na kraju treba da imamo konacnu matricu H koja nije 3D matrica, vec 2D
    # Isto tako i za R
    # I na kraju ide ovaj 5. korak - ali je moj mozak vec sprzen do kraja, pa nemam snage da ovo dalje kucam
    # Inace bilo bi pametno proci i provjeriti da li nam se dimenzije svih matrica slazu
    # Mislim kad se radi mnozenje, sabiranje i to
    minima = d.min(axis=1) # Minimumi svake vrste
    
    print 'Minimumi u matrici d: ', minima
    map_index = d.argmin(axis=1) # Sama vrijednost ove matrice predstavlja ono gornje j neko, a indeks u listi samo i odozgo
    print 'Indeksi tih minimuma: ', map_index
    g = 0.5 # Rekao je da probamo sa 0.5, 1 i 5, pa da vidimo kako se ponasa
    measurement_index = np.nonzero(minima<g**2)
    print 'Measurement index: ', measurement_index
    map_index = map_index[measurement_index[0]] # Indeksi po vrstama matrice d, koji su valjani
    print 'Measurement index[0]',measurement_index[0] 
    print 'Novo map_index: ', map_index
    print 'Nmapentries: ', nMapEntries
    print 'ovo nesto: ', map_index + measurement_index[0]*nMapEntries
    
    
    
    v_ino = v_ino[:, map_index[0] + measurement_index[0][0]*nMapEntries]
   
    print 'V na kraju: ', np.shape(v_ino)
    H = H[:,:,map_index[0]] # Imamo shit sa ovim indeksima ovdje
    print 'H na kraju: ', np.shape(H)

	
    # 5. korak: Estimacija
    print '========4.korak - ========'
    # S je krajnja kovarijansa krajnje inovacije
    #S = H * P_priori * H.T + R
    S = np.linalg.multi_dot([H,P_priori,H.T])+R
    # Kalmanovo pojacanje
    #K = P_priori * H.T*np.linalg.inv(S)
    K = np.linalg.multi_dot([P_priori,H.T,np.linalg.inv(S)])
    print 'K', np.shape(K)
    P_posteriori = np.dot((np.eye(3) - np.dot(K,H)),P_priori)
    x_posteriori = x_priori + np.dot(K,v_ino)
	
    x = x_posteriori
    P =P_posteriori
    print 'xNOVO: ', x
    print 'P novo: ', P
    a = Float32MultiArray(data=x)
    # I onda na kraju je potrebno publishovati ovo x_posteriori, i pored toga ono Z, pa da poredimo to mjerenje
    pub.publish(a)

def timer_callback(event):
    # print 'Timer called at ' + str(event.current_real) # ovo je ako hocemo da testiramo da li radi tajmer 
    # pozivamo main funkciju, ovo bi sad trebalo da se poziva na svake 0.2 sekunde
    global rho_temp, phi_temp
    glavna(rho_temp, phi_temp) # mozemo i main da podesimo odmah kao callback funkciju ali ajdeee... 
    print 'Srki PROBUDI SE'
def callback_odom(odom_data):
    # Najbolje je sa odoma da citamo ostvarene linearne i ugaone brzine, i preko njih nadjemo delta_s_r i delta_s_l
    global v
    v = 0
    print 'Brzina v u odomu', v
    global w
    global vl 
    global vr
	
    l = 0.16/2 # Polovina rastojanja izmedju tockova (u specifikacijama)
    print l
	
    # Ne znam kako da ucitam brzine iz Twist-a:
    # Ovde je objasnjeno  nesto:
    # https://stackoverflow.com/questions/50976281/what-are-the-x-y-and-z-mean-in-ros-geometry-msgs-twist <- potencijalno zanemariti postojanje ovog linka
    # http://wiki.ros.org/navigation/Tutorials/RobotSetup/Odom <-- ovaj link je relevantniji
    vx = odom_data.twist.twist.linear.x
	
    # vy = odom_data.twist.twist.linear.y
    # vz = 0
	
    wz = odom_data.twist.twist.angular.z
    print 'Brzina vx u clbodom: ', vx
    print 'Brzina wz u clbodom: ', wz
	
    v = vx;
    w = wz;
	
    # ili:
    # v = math.sqrt(vx^2 + vy^2 + vz^2); # ovo sve treba proveriti ofc
	
    vl = v - w*l
    vr = v + w*l
		
def callback_laser(laser_data):
    # Ovdje ide kod koji mi samo ucitava podatke sa Lasera, znaci samo da ucitam ono rho_temo i phi_temo, i poslije
    # ih obradjujem u glavnoj funkciji, tako da ugl sve one funkcije treba kopirati iz prethodnog zadatka, ali mislim
    # da se to moze odraditi naknadno, da ne bismo sad imali 1000 linija koda bezveze, znamo da taj dio koda samo treba da izbaci krajnje podatke
    # u oblikeu M = np.vstack((alpha, r)), gdje prvu vrstu ove matrice cine sve alphe, a drugu vrstu cine svi r-ovi
    global rho_temp
    global phi_temp
    global world_map
    global first
    rho_temp = laser_data.ranges
    print 'Inkrement', laser_data.angle_increment
    phi_temp = []
    phi_temp = list(np.arange(laser_data.angle_min,6.29,laser_data.angle_increment))
  
    print 'Rho_temp u CLBKLASER: ', len(rho_temp)
    print 'Phi_temo u CLBKLASER: ', len(phi_temp)
	
    if (first == 1): 
	# Ovdje je potrebno kopirati onaj kod iz prethodnog zadatka da bi se dobili podaci o mapi
	#!!!!!!
        alpha, r = SAM(rho_temp, phi_temp)
	world_map = np.vstack((alpha,r)) 
	first = 0
	
def listener():
    # Inicijalizujem node
    rospy.init_node('domaci4', anonymous = True)
    # Definisem prvi subscriber koji ce obradjivati podatke u prvom callbacku
    rospy.Subscriber('/odom', Odometry, callback_odom) 
    # Definisem drugi subscriber koji ce obradjivat podatke u drugom callbacku
    rospy.Subscriber('scan', LaserScan, callback_laser)
    rospy.Timer(rospy.Duration(4), timer_callback, oneshot = False)
    # Dodavanje tajmera !!!
    # kako se koristi ali je kod u c++: http://answers.ros.org/question/199727/how-to-use-timer/
    # kako se poziva: http://wiki.ros.org/rospy/Overview/Time
    # oneshot = false znaci da ce tajmer periodicno da se izvrsava
    # callback funkcija se poziva na svakih rospy.Duration sekundi
    # proveriti da li moze ovo u listeneru ili treba da bude negde drugde !!!
    
    # Pitati ga za ovaj tajmer
    rospy.spin()

if __name__=='__main__':
    global rho_temp
    global phi_temp
    global v
    v = 0
    global w
    w = 0
    global x
    global vl
    vl = 0
    global vr
    vr = 0
    global first
    first = 1
    print 'Flag: ', first
    global world_map
    x = [0,0,0] # Pretpostavljeno pocetno stanje robota je 0,0 i valjda isto rotacija je 0
    print 'Pocetno stanje: ', x
    global P	
    P = np.array([[0.0001,0,0],[0,0.0001,0],[0,0,0.0001]])
    print 'P0 iznosi: ', P
    
    listener() 
