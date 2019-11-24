#! /usr/bin/env python

import rospy
from nav_msgs.msg import Odometry
from tf.transformations import euler_from_quaternion
from geometry_msgs.msg import Point, Twist
import math

# Funkcija za normalizaciju ugla
def normalizeAngle(angle):
    angle1 = (angle+math.pi)%(2*math.pi)-math.pi
    return angle1

# Kontroler
def controllerOutput(x, y, theta, x_g, y_g, theta_g, backward, krho, kalpha, kbeta, constSpeed):
    d_x = x_g-x
    d_y = y_g-y

    # Racunanje kontrolnih velicina
    rho = math.sqrt(d_x**2+d_y**2)
    lambd = math.atan2(d_y, d_x)
    alpha = lambd-theta
    alpha = normalizeAngle(alpha)

    if backward: # ukoliko je dozvoljeno da se robot krece i unazad
        if (abs(alpha)<=(math.pi/2)):
            beta = theta_g-lambd
            Krho = krho
        else:
            alpha = lambd-theta-math.pi
            alpha = normalizeAngle(alpha)
            beta = theta_g-lambd-math.pi
            Krho = -krho
    else:
        beta = theta_g-lambd
        Krho = krho

    beta = normalizeAngle(beta)
    v = Krho*rho
    omega = kalpha*alpha+kbeta*beta

    absVel = abs(v)
    if (absVel>10**(-6)):
    	v = v/absVel*constSpeed
    	omega = omega/absVel*constSpeed

    # Na kraju se vraca v i omega
    lin_x = v*math.cos(theta)
    lin_y = v*math.sin(theta)

    return lin_x, lin_y, omega

def callback(msg):
    global x
    global y
    global theta
    global EndCond

    # Inicijalizacija strukture koju se salje preko teme cmd_vel
    speed = Twist()

    # Stvaranje objavljivaca
    pub = rospy.Publisher("/cmd_vel", Twist, queue_size = 1)

    if (not(EndCond)):
	    # Ucitavanje trenutnih koordinata i orijentacije robota
            x = msg.pose.pose.position.x
            y = msg.pose.pose.position.y
            rot_q = msg.pose.pose.orientation
            (roll, pitch, theta) = euler_from_quaternion([rot_q.x, rot_q.y, rot_q.z, rot_q.w])

            speed.linear.x, speed.linear.y, speed.angular.z = controllerOutput(x, y, theta, x_g, y_g, theta_g, backward, krho, kalpha, kbeta, constSpeed)

            d_theta = abs(normalizeAngle(theta-theta_g))
            rho = math.sqrt((x_g-x)**2+(y_g-y)**2)
            EndCond = (rho<dist_threshold) and (d_theta<angle_threshold)

            pub.publish(speed)

def listener():
	rospy.init_node('controller_4', anonymous=True)
	rospy.Subscriber("/odom", Odometry, callback)
	rospy.spin()

if __name__=='__main__':
    try:
	# Inicijalizacija
        x = 0.0
        y = 0.0
        theta = 0.0
	EndCond = 0

	# Ucitavanje koordinata cilja i ugla robota
    	x_g = float(input("Unesi koordinatu x [m]:"))
    	y_g = float(input("Unesi koordinatu y [m]:"))
    	theta_g = float(input("Unesi ugao theta_g [rad]:"))
    	backward = int(input("Dozvoli kretanje unazad (0/1):"))
        constSpeed = float(input("Unesi vrijednost brzine (0-1 m/s):"))
	
	# Ucitavanje vrijednosti pragova za zaustavljanje robota
    	dist_threshold = 0.25
    	angle_threshold = 0.1 # oko 5 stepeni

    	# Ucitavanje parametara kontrolera
    	krho = 0.5
    	kalpha = 1.5
    	kbeta = -0.6

        listener()
    except rospy.ROSInterruptException:
		pass
