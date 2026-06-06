from tkinter import *
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.figure
from mpl_toolkits.mplot3d import Axes3D
import tkinter.messagebox as tmsg
import numpy as np

#####################################################################################################
########################################################################################################
class GUI(Tk):
    def __init__(self):
        super().__init__()
        #this fixes the size
        self.geometry("1100x640")
        self.maxsize(1100,640)
        self.minsize(1100,640)
        #title of our gui
        self.title("Observer and Celestial Coordinates")

###################################################################################################################
    def plot(self, *args):
        latitude = self.latitude.get()*(np.pi/180)
        time = self.time.get()
        day = self.day.get()

        dec = self.dec.get()*(np.pi/180)
        alpha = self.alpha.get()*15.0*np.pi/180

        # Create a meshgrid for the sphere
        phi_osphere, theta_osphere = np.mgrid[0.0:2.0*np.pi:100j, 0.0:np.pi:50j]
        x_osphere = np.sin(theta_osphere) * np.cos(phi_osphere)
        y_osphere = np.sin(theta_osphere) * np.sin(phi_osphere)
        z_osphere = np.cos(theta_osphere)

        #get the sun's parameters based on a straight line fit. go from 0 to 24 in 360 days
        alpha_sun = 0.0 + ((24.0 - 0.0)/(360.0 - 1.0))*(day - 1.0)
        alpha_sun = alpha_sun*15.0*np.pi/180.0
        if 1<= day <=90: #spring to summer (weather convention for the northern hemisphere. Sorry Aussies!)
            dec_sun = 0.0 + ((23.5 - 0.0)/(90.0 - 1.0))*(day - 1.0)
        elif 91<= day <= 270:#summer to winter
            dec_sun = 23.5 + ((-23.5 - 23.5)/(270.0 - 90.0))*(day - 90.0)
        elif 271<=day<=360: #winter to spring
            dec_sun = -23.5 + ((0.0 - (-23.5))/(360.0 - 270.0))*(day - 270.0)
        dec_sun = dec_sun*np.pi/180.0
        HA_sun = time - 12.0 #number of hours past local noon
        HA_sun = HA_sun*15.0*np.pi/180.0
        phi_sun = - HA_sun

        #for the test object
        HA = HA_sun + alpha_sun - alpha
        phi = - HA

        #radius of sphere
        R = 1.0
        
        #unrotated 
        #the circle of the sun
        zsun_day_old = np.repeat(R*np.sin(dec_sun),100)
        xsun_day_old = R*np.cos(dec_sun)*np.cos(np.linspace(0,2.0*np.pi,100))
        ysun_day_old = R*np.cos(dec_sun)*np.sin(np.linspace(0,2.0*np.pi,100))
        #the point of the sun
        zsun_old = R*np.sin(dec_sun)
        xsun_old = R*np.cos(dec_sun)*np.cos(phi_sun)
        ysun_old = R*np.cos(dec_sun)*np.sin(phi_sun)
        #the line
        zpole_old = np.linspace(-1.0,1.0,100)
        xpole_old = np.repeat(0,100)
        ypole_old = np.repeat(0,100)
        #the circle of test object
        z_day_old = np.repeat(R*np.sin(dec),100)
        x_day_old = R*np.cos(dec)*np.cos(np.linspace(0,2.0*np.pi,100))
        y_day_old = R*np.cos(dec)*np.sin(np.linspace(0,2.0*np.pi,100))
        #the point of the test object
        z_old = R*np.sin(dec)
        x_old = R*np.cos(dec)*np.cos(phi)
        y_old = R*np.cos(dec)*np.sin(phi)

        theta_sun = -(np.pi/2.0 - latitude) #rotation matrix angle
        theta = - (np.pi/2.0 - latitude)
        #rotated coordinates
        xsun_day = np.cos(theta_sun)*xsun_day_old + np.sin(theta_sun)*zsun_day_old
        ysun_day = ysun_day_old
        zsun_day = -np.sin(theta_sun)*xsun_day_old + np.cos(theta_sun)*zsun_day_old

        xsun = np.cos(theta_sun)*xsun_old + np.sin(theta_sun)*zsun_old
        ysun = ysun_old
        zsun = -np.sin(theta_sun)*xsun_old + np.cos(theta_sun)*zsun_old

        xpole = np.cos(theta)*xpole_old + np.sin(theta)*zpole_old
        ypole = ypole_old
        zpole = -np.sin(theta)*xpole_old + np.cos(theta)*zpole_old

        x_day = np.cos(theta)*x_day_old + np.sin(theta)*z_day_old
        y_day = y_day_old
        z_day = -np.sin(theta)*x_day_old + np.cos(theta)*z_day_old

        x = np.cos(theta)*x_old + np.sin(theta)*z_old
        y = y_old
        z = -np.sin(theta)*x_old + np.cos(theta)*z_old

        #the observer frame plot
        self.ax.cla()
        self.ax.set_title('Observer Frame')
        self.ax.set_axis_off()
        self.ax.scatter(xsun,ysun,zsun, c='red', s= 50, label='Sun')
        self.ax.plot(xsun_day,ysun_day,zsun_day,c='red', linestyle='dashed')
        self.ax.scatter(x,y,z, c='blue', s= 50, label=' Test Oject')
        self.ax.plot(x_day,y_day,z_day,c='blue', linestyle='dashed')
        self.ax.plot(xpole,ypole,zpole,c='green', linestyle='dashed', label='Earth Axis')

        self.ax.plot(1.0*np.cos(np.linspace(0,2.0*np.pi,100)), 1.0*np.sin(np.linspace(0,2.0*np.pi,100)),np.zeros(100), c='black')
        self.ax.plot(np.linspace(-1,1,100), np.zeros(100), np.zeros(100), c='black', linestyle='dashed')
        self.ax.plot(np.zeros(100), np.linspace(-1,1,100), np.zeros(100), c='black', linestyle='dashed')
        self.ax.plot(np.zeros(100),np.zeros(100),np.linspace(-1,1,100),c='black', linestyle='dashed')
        self.ax.plot_surface(x_osphere, y_osphere, z_osphere, color='b', alpha=0.2)
        
        self.ax.text(-1.3,0,0, "N", color='green', fontsize=12, ha='center', va='center')
        self.ax.text(0,0,1.3, "zenith", color='green', fontsize=12, ha='center', va='center')
        self.ax.text(0,0,-1.3, "nadir", color='green', fontsize=12, ha='center', va='center')
        self.ax.view_init(5, 250)

        self.ax.legend(loc=3)
        
        #the celestial coordinates plot
        #straight line to be at phi=0 at day 90/summer solstice
        phi_earth = 0.0 + ((3.0*np.pi/2.0 - (-np.pi/2.0))/(360.0-1.0))*(day-90.0)
        phi_object = alpha + np.pi/2.0 #since 0h is at phi = 90 degrees
        self.ax2.cla()
        self.ax2.set_title('Celestial Frame')
        self.ax2.set_axis_off()
        self.ax2.scatter(0,0,0, c='red', s= 100, label='Sun')
        self.ax2.scatter(1.0*np.cos(phi_earth),1.0*np.sin(phi_earth),0, c='green', s= 50, label='Earth')
        self.ax2.scatter(1.0*np.cos(dec)*np.cos(phi_object), 1.0*np.cos(dec)*np.sin(phi_object),1.0*np.sin(dec), c='blue', label='Test Object')
        self.ax2.plot(1.0*np.cos(phi_earth) + np.linspace(-0.25*np.tan(23.5*np.pi/180.0),0.25*np.tan(23.5*np.pi/180.0),100), 1.0*np.sin(phi_earth)+ np.zeros(100), 0 - np.linspace(-0.25,0.25,100), c='green',label='Earth Axis')
        self.ax2.plot(1.0*np.cos(phi_earth) + np.zeros(100), 1.0*np.sin(phi_earth)+ np.linspace(-0.5,0.5,100), np.zeros(100), c='purple',linewidth=1, label='First Point of Aries')
        self.ax2.plot(1.0*np.cos(np.linspace(0,2.0*np.pi,100)), 1.0*np.sin(np.linspace(0,2.0*np.pi,100)),np.zeros(100), c='black')
        self.ax2.plot(np.linspace(-1,1,100), np.zeros(100), np.zeros(100), c='black', linestyle='dashed')
        self.ax2.plot(np.zeros(100), np.linspace(-1,1,100), np.zeros(100), c='black', linestyle='dashed')
        self.ax2.plot_surface(x_osphere, y_osphere, z_osphere, color='b', alpha=0.2)
        self.ax2.set_zlim([-1,1])
        self.ax2.view_init(20, 250)
        self.ax2.text(1.3,0,0, "18 h", color='green', fontsize=12, ha='center', va='center')
        self.ax2.text(-1.3,0,0, "6 h", color='green', fontsize=12, ha='center', va='center')
        self.ax2.text(0,1.3,0, "0 h", color='green', fontsize=12, ha='center', va='center')
        
        self.ax2.legend(loc=3)



        self.canvas.draw()
#############################################################################################################################################
    def make_title(self):
        f1 = Frame(self, bg='red', borderwidth=10, relief=GROOVE)
        f1.pack(side=TOP, fill=X, pady=10)
        Label(f1, bg = "red", fg = "black", font = ("calibiri", 20, "bold"),  text = "Observer and Celestial Coordinates").pack()
######################################################################################################################################
    def make_selection_frame(self):
        f2 = Frame(self, bg = "lightblue", borderwidth=10 , relief = RAISED)
        f2.pack(fill=X)

        Label(f2, bg = "lightblue", fg = "black", font = ("calibiri", 15, "bold"), text = "Choose Location, Time and Day. Set Test Object in Frame Below").pack()
        
        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Latitude (Deg)").pack(side = LEFT, padx=10)
        self.latitude = Scale(f2, bg ='lightblue', fg='black',from_=-90, to=90, orient=HORIZONTAL, resolution=1, length = 200)
        self.latitude.set(0)
        self.latitude.pack(side = LEFT, padx=10)

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text=" Solar Time (h)").pack(side = LEFT, padx=10)
        self.time = Scale(f2, bg ='lightblue', fg='black',from_=0, to=24, orient=HORIZONTAL, resolution=0.5, length=200)
        self.time.set(12)
        self.time.pack(side = LEFT, padx=10)

        Label(f2, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Day").pack(side = LEFT, padx=10)
        self.day = Scale(f2, bg ='lightblue', fg='black',from_=1, to=360, orient=HORIZONTAL, resolution=1, length=400)
        self.day.set(0)
        self.day.pack(side = LEFT, padx=10)
        
        self.latitude.config(command = self.plot)
        self.time.config(command = self.plot)
        self.day.config(command = self.plot)

############################################################################################################################################
    def make_plotframe(self):
        f3 = Frame(self, bg = "lightblue", borderwidth = 10, relief = RAISED)
        f3.pack(fill=X)

        self.fig = matplotlib.figure.Figure(edgecolor='black', facecolor='red',linewidth=7)
        self.fig.set_size_inches(8,4, forward=True)
        self.ax = self.fig.add_subplot(1,2,1, projection='3d')
        self.ax2 = self.fig.add_subplot(1,2,2, projection='3d')
        #self.ax2.set_facecolor('black')
        
        self.fig.subplots_adjust(right=1)
        self.fig.subplots_adjust(left=0)
        
        Label(f3, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Declination (Deg)").pack(side = LEFT, padx=10)
        self.dec = Scale(f3, bg ='lightblue', fg='black',from_=-90, to=90, orient=VERTICAL, resolution=1, length=200)
        self.dec.set(21)
        self.dec.pack(side=LEFT, padx=10)

        Label(f3, bg = "lightblue", fg ="black", font = ("calibiri", 10, "bold"), text="Right Ascension (h)").pack(side = RIGHT, padx=10)
        self.alpha = Scale(f3, bg ='lightblue', fg='black',from_=0, to=24, orient=VERTICAL, resolution=0.5, length=200)
        self.alpha.set(6)
        self.alpha.pack(side = RIGHT, padx=10)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=f3)  
        self.canvas.get_tk_widget().pack( pady=5)

        # Create Toolbar
        toolbar = NavigationToolbar2Tk(self.canvas, f3, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(fill=X) 

        self.dec.config(command = self.plot)
        self.alpha.config(command = self.plot)

#####################################################################################################
        
#The main body of the program, which just calls all the functions from before

#instance of out GUI class
observer_celestial = GUI()

#methods that give life to the GUI
observer_celestial.make_title()
observer_celestial.make_selection_frame()
observer_celestial.make_plotframe()
observer_celestial.plot()

#the infinite GUI loop
observer_celestial.mainloop()