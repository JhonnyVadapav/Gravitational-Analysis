import warnings
warnings.filterwarnings('ignore')
try:
    print('\t\t\t ----------------------')
    print('\t\t\t GRAVITATIONAL ANALYSIS')
    print('\t\t\t ----------------------')
    print('''Efforts by-
    Harshit Arora
    Raghav Sakhuja
    Siddharth Verma''')
    print()
    print()
    print()
    print('DISCLAIMER!!!')
    print('''Please make sure to install the following modules to run the given codes-
    # matplotlib
    # astropy
    # astroquery
    # plotly
    ''')
    while True:
        print()
        print('MENU-')
        print('-'*67)
        x=' '
        print('| #1  | Simulation ',x*45,'|')
        print('| #2  | Graph for |g|',x*43,'|')
        print('| #3  | Database',x*48,'|')
        print('| #4  | Variation in gravitational field due to rotation',x*8,'|')
        print('| #5  | Gravitational potential energy of a planet',x*14,'|')
        print('| #6  | Max possible structural height on a planet',x*14,'|')
        print('| #7  | Orbital velocity at a height h above surface of a planet','','|')
        print('| #8  | Escape velocity from the surface of planet',x*14,'|')
        print('| #9  | Time period of a satellite revolving around a planet',x*4,'|')
        print('| #10 | Altitude of Kármán Line',x*33,'|')
        print('| #11 | Gravitational redshift',x*34,'|')
        print('| #12 | Exit',x*52,'|')
        print('-'*67)
        print()
        ch=int(input('Enter the # of your choice:'))
        lch=[1,2,3,4,5,6,7,8,9,10,11,12]
        if ch not in lch:
            print()
            print('Please enter a valid choice')
            print()
            continue
        elif ch==1:
            print('\t\t\t ----------')
            print('\t\t\t Simulation')
            print('\t\t\t ----------')
            print()
            while True:
                print('''
    |#1| Begin the simulation
    |#2| To the menu''')
                j=int(input('Enter your choice:'))
                if j==1:
                    print()
                    print('''This simulation will show you the movement
    of the first 4 planets around the sun.''')
                    print('''
    The simulation is in real time and starts at 2018-01-01
    and stops at 2020-01-04.
    Only 4 planets are shown because the distancce of Jupiter and the other
    outer planets is large enough to make the maneuverability in the
    simulation difficult.''')
                    print('#The simulation is to scale except for the sun#')
                    print()
                    print('Shift to full screen for a smoother simulation.')
                    print()
                    start=input('Input start to begin:')
                    import numpy as np
                    import matplotlib.pyplot as plt
                    import matplotlib.animation as animation
                    from astropy.time import Time
                    from astroquery.jplhorizons import Horizons

                    sim_start_date = "2018-01-01"     # simulating a solar system starting from this date
                    sim_duration = 2 * 365                # (int) simulation duration in days
                    m_earth = 5.9722e24 / 1.98847e30  # Mass of Earth relative to mass of the sun
                    m_moon = 7.3477e22 / 1.98847e30

                    class Object:                   # define the objects: the Sun, Earth, Mercury, etc
                        def __init__(self, name, rad, color, r, v):
                            self.name = name
                            self.r    = np.array(r, dtype=np.float)
                            self.v    = np.array(v, dtype=np.float)
                            self.xs = []
                            self.ys = []
                            self.plot = ax.scatter(r[0], r[1], color=color, s=rad**2, edgecolors=None, zorder=10)
                            self.line, = ax.plot([], [], color=color, linewidth=1.4)

                    class SolarSystem:
                        def __init__(self, thesun):
                            self.thesun = thesun
                            self.planets = []
                            self.time = None
                            self.timestamp = ax.text(.03, .94, 'Date: ', color='w', transform=ax.transAxes, fontsize='x-large')
                        def add_planet(self, planet):
                            self.planets.append(planet)
                        def evolve(self):           # evolve the trajectories
                            dt = 1.0
                            self.time += dt
                            plots = []
                            lines = []
                            for p in self.planets:
                                p.r += p.v * dt
                                acc = -2.959e-4 * p.r / np.sum(p.r**2)**(3./2)  # in units of AU/day^2
                                p.v += acc * dt
                                p.xs.append(p.r[0])
                                p.ys.append(p.r[1])
                                p.plot.set_offsets(p.r[:2])
                                p.line.set_xdata(p.xs)
                                p.line.set_ydata(p.ys)
                                plots.append(p.plot)
                                lines.append(p.line)
                            self.timestamp.set_text('Date: ' + Time(self.time, format='jd', out_subfmt='date').iso)
                            return plots + lines + [self.timestamp]

                    plt.style.use('dark_background')
                    fig = plt.figure(figsize=[6, 6])
                    ax = plt.axes([0., 0., 1., 1.], xlim=(-1.8, 1.8), ylim=(-1.8, 1.8))
                    ax.set_aspect('equal')
                    ax.axis('off')
                    ss = SolarSystem(Object("Sun", 40, 'yellow', [0, 0, 0], [0, 0, 0]))
                    ss.time = Time(sim_start_date).jd
                    colors = ['gray', 'orange', 'blue', 'chocolate','white']
                    sizes = [0.069, 0.17, 0.2, 0.097]
                    names = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter']
                    texty = [0.5, 0.75, 1.075, 1.5]
                    for i, nasaid in enumerate([1, 2, 3, 4]):  # The 1st, 2nd, 3rd, 4th planet in solar system
                        obj = Horizons(id=nasaid, location="@sun", epochs=ss.time, id_type='id').vectors()
                        ss.add_planet(Object(nasaid, 20 * sizes[i], colors[i], 
                                             [np.double(obj[xi]) for xi in ['x', 'y', 'z']], 
                                             [np.double(obj[vxi]) for vxi in ['vx', 'vy', 'vz']]))
                        ax.text(0, - (texty[i] + 0.1), names[i], color=colors[i], zorder=1000, ha='center', fontsize='large')
                    def animate(i):
                        return ss.evolve()
                    ani = animation.FuncAnimation(fig, animate, repeat=False, frames=sim_duration, blit=True, interval=20,)
                    plt.show()
                elif j==2:
                    break
        elif ch==2:
            print('\t\t\t -------------')
            print('\t\t\t Graph for |g|')
            print('\t\t\t -------------')
            print('''      This program shows the variation in the strength of the
          gravitational acceleration i.e. |g| with increasing
          distance from the center of the planetary body by
          plotting a graph between |g| and the distance r.''')
            print()
        
            '''Matplotlib is a plotting library for the Python programming language and
                its numerical mathematics extension NumPy. Using matplotlib as it creates a
                figure, creates a plotting area in a figure, plots some lines in a plotting
                area, decorates the plot with labels, etc.'''
            print('''|g| for two bodies/objects can be calculated using the following formula-
                \t        |g|= G(M1)(M2)/R^2
            Where,
            G  ---> Universal Gravitational Constant
            M1 ---> Mass of first object
            M2 ---> Mass of second object
            R  ---> Distance between the objects''')
            print('\n')
            import matplotlib.pyplot as plt
            l=[]
            G=6.67*(10**-11)
            while True:
                radius_input=int(input('Enter the the significant figures in radius of the planet in m  :'))
                exp1=int(input('Enter the exponent of 10 used in the radius  :'))
                power1=10**(exp1)
                R=radius_input*power1
                print('Preview of the radius:')
                print(R,'m',sep='')
                j=input('''
    If the preview is the desired radius then input yes
    otherwise input no:''')
                if j=='no':
                    continue
                elif j=='yes':
                    break
            print()
            y=input('''
    # If mass of the planet is known, enter M
    # If average density is known, enter p:''')
            print()
            if y=='M':
                print()
                print('''# Please enter the mass up to second decimal place''')
                while True:
                    print()
                    m=float(input('Enter the mass:'))
                    e=int(input('Enter the power of 10:'))
                    print('Preview of the mass:')
                    M=m*(10**e)
                    print(M,'kg',sep='')
                    ch=input('If the preview is the desired mass then input yes otherwise input no:')
                    if ch=='yes':
                        M=m*(10**e)
                        break
                    elif ch=='no':
                        continue
                g=G*M/R**2
            elif y=='p':
                p=float(input('Enter average density of the planet(in kg/m^3):'))
                g=(4*G*R*p*3.14)/3
            print()
            print('The value of g on the surface=',g)
            '''This user defined function calculates the value of |g| when the value
            the distance from the centre of the planet is more than the radius of the planet'''
            def gR(Ro):
                G=6.67*(10**-11)
                if y=='M':
                    g_out=G*M/Ro**2
                elif y=='p':
                    m=(4*3.14*p*R**3)/3
                    g_out=G*m/Ro**2
                return g_out
            '''This user defined function calculates the value of |g| when the value
               the distance from the centre of the planet is less than the radius of the planet'''
            def gr(r):  #a point inside the surface 
                g_in=(g*r)/R
                return g_in
            print()
            r=int(input('Initial distance from the centre(in metres):'))
            print()
            x=int(input('Final distance from the centre(in metres):'))
            for i in range(r,x+1,100000):
                if i<R:
                    l.append(gr(i))
                elif i>=R:
                    l.append(gR(i))
            '''the values of |g| are stored as a list which are then used to calculate the graph 
               for the graph'''
            for k in l:
                print(k)
            print('The above are the values of |g| at various step values entered before.')
            print()
            print('''If a '|g| vs r' graph is plotted then it will show a linear relationship
    between |g| and r till value of r is less than the radius of the planet. When the
    value of r becomes larger than the radius of the planet, the graph follows the
    relation y=1/(x^2) hence a curved line.''')
            print('\n')
            u=input('Do you want to see the graph(yes/no):')
            if u=='yes':
                #with the matploylib function we make the graph for values of |g| 
                plt.plot(l)
                #r/10^5 values on the x-axis
                plt.xlabel('r/10^5 --->',fontsize=20,color='black') 
                #|g| values on the y-axis
                plt.ylabel('|g| --->',fontsize=20,color='black') 
                plt.title('graphical representation of g at different intervals')
                ax=plt.axes()
                ax.set_facecolor("black")
                plt.savefig('myplot.png') 
                plt.show()
                print()
                print('#1 To the menu')
                print('#2 Terminate the program')
                j=int(input('Enter your choice(1/2):'))
                if j==1:
                    continue
                if j==2:
                    print('Program terminated')
                    break
            elif u=='no':
                print()
                print('#1 To the menu')
                print('#2 Terminate the program')
                j=int(input('Enter your choice(1/2):'))
                if j==1:
                    continue
                if j==2:
                    print('Program terminated')
                    break
            continue
            print()
            print('#1 To the menu')
            print('#2 Terminate the program')
            j=int(input('Enter your choice(1/2):'))
            if j==1:
                continue
            if j==2:
                print('Program terminated')
                break
            continue
        elif ch==3:
            print('\t\t\t --------')
            print('\t\t\t Database')
            print('\t\t\t --------')
            print("Which databases would you like to view")
            print("(please make sure you have run the codes 4 , 5 , 9 before viewing the databases 2 , 3 , 4 )")
            print("#1 : Planets")
            print("#2 : Time period of a satellite")
            print("#3 : g-field due to rotation of Earth")
            print("#4 : GPE of a planet")
        
            
            c=int(input("Enter your choice:"))
            if c==1:
                import plotly.graph_objects as go
                headerColor = 'grey'
                rowEvenColor = 'lightgrey'
                rowOddColor = 'white'
                fig = go.Figure(data=[go.Table(
                    header=dict(
                        values=['<b>PLANETS</b>','<b>MASS(in kg)</b>','<b>RADIUS(in m)</b>','<b>NUMBER OF MOONS</b>','<b>g value</b>','<b>ECCENTRICITY</b>','<b>DISTANCE FROM SUN(AU)</b>'],
                        line_color='darkslategray',
                        fill_color=headerColor,
                        align=['left','center'],
                        font=dict(color='white', size=12)
                    ),
                    cells=dict(
                        values=[
                            ['MERCURY', 'VENUS', 'EARTH', 'MARS','JUPITER','SATURN','URANUS','NEPTUE'],
                            ['0.30x10^24','4.87x10^24','5.97x10^24','0.642x10^24','1898x10^24','568x10^24','86.8x10^24','102x10^24'],
                            ['2.43x10^6','6.06x10^6','6.37x10^6','3.37x10^6','6.99x10^7','5.85x10^7','2.33x10^7','2.21x10^7'],
                            ['0','0','1','2','79','62','27','14'],
                            ['3.61','8.83','9.81','3.75','26.0','11.2','10.5','13.3'],
                            ['0.206','0.007','0.017','0.093','0.048','0.056','0.047','0.009'],
                            ['0.39','0.723','1','1.524','5.203','9.539','19.18','30.06']],
                        line_color='darkslategray',
                        # 2-D list of colors for alternating rows
                        fill_color = [[rowOddColor,rowEvenColor,rowOddColor, rowEvenColor,rowOddColor]*5],
                        align = ['left', 'center'],
                        font = dict(color = 'darkslategray', size = 11)
                        ))
                ])
                fig.show()
                
            elif c==2 : 
                print("[radius of planet , Mass of planet , time period]")
                f=open("gravity1.text",'r') 
                print(f.read())
                f.close()
        
            elif c==3 :
               print("[ Radius , Mass , g-pole , g-equator ]")
               f=open("gfield.text",'r') 
               print(f.read())
               f.close()
                    
            elif c==4 :
              print("[ Radius of planet , Mass , GPE ]")
              f=open("potential1.text",'r') 
              print(f.read())
              f.close()
            print()
            print('#1 To the menu')
            print('#2 Terminate the program')
            j=int(input('Enter your choice(1/2):'))
            if j==1:
                continue
            if j==2:
                print('Program terminated')
                break
            continue
            

        elif ch==4:
            print('\t\t\t ------------------------------------------------')
            print('\t\t\t Variation in gravitational field due to rotation')
            print('\t\t\t ------------------------------------------------')
            file=open("gfield.text",'a')
            record=[]
            r=[]
            
            print('''
     Let us consider the earth to be a spherical ball of mass ‘M’ and radius ‘R’.
     An object of mass ‘m’ is at point P at latitude φ, when the earth is not rotating the weight of the object is mg.
     But earth is rotating with angular velocity ω. So, the object is moving in a circular path of radius ‘r’ as shown in the figure.
     The object experiences centrifugal force,

                  F=mw^2r=mw^2rcosQ
                  
    The object is being acted by two forces ‘mg’ and ‘F‘. The resultant of these two forces gives the apparent weight of the object (m’g).
    on solving this we get
                  g'=sqrt(g^2+(w^2r)^2-2grw^2sinQ

    at equater Q=90*
                  thus g'=g-w**2r
    at poles Q=0
                  thus g'=g
    therefore , at poles we feel the maximum acceleration due to gravity
            '''
            )
            print()
            import math
            import matplotlib.pyplot as plt
            G=6.67*(10**-11)
            glist=[]
            def g_rotation(angle,speed,g,r):
                g_due_rot=(g*g-2*(speed**2)*(math.sin(angle))*r*g+(speed*speed*r)**2)**0.5
                return g_due_rot
            #function to find mass of planet using density
            def mass(radius,rho):
                m_o_planet=rho*(4*3.14*(radius**3))/3
                return m_o_planet
            #we input the values
            while True:
                radius_input=int(input('enter the the significant figures in radius of the planet in m  :'))
                exp1=int(input('enter the exponent of 10 used in the radius  :'))
                power1=10**(exp1)
                R_of_planet=radius_input*power1
                r.append(R_of_planet)
                print('the required radius is:',R_of_planet)
                j=input('''If the preview is the desired radius then input yes otherwise input no:''')
                if j=='no':
                    continue
                elif j=='yes':
                    break
            print('if you know the mass of the planet , enter M')
            print('if you know the density of the planet enter p')
            option=input('enter your choice  :')
            if option=='M':
                #we input significant of mass and the power separately to reduce the work of the user
                while True:
                    mass1=int(input('enter the significant figures of mass of the planet  :'))
                    exp=int(input('enter the exponent of 10 used in the mass  :'))
                    power=10**(exp)
                    Mass=mass1*power
                    print('the required mass is:',Mass)
                    j=input('''If the preview is the desired mass then input yes otherwise input no:''')
                    if j=='no':
                        continue
                    elif j=='yes':
                        break
            elif option=='p':
                while True:
                    print('''We find the mass of the planet using the Density of the planet using the formula
                    (4*3.14*(Radius**3)/3)*rho''')
                    rho_planet=int(input('enter the density of the planet in kg/m3'))
                    Mass=mass(R_of_planet,rho_planet)
                    print('the required mass is:',Mass)
                    j=input('''If the preview is the desired mass then input yes otherwise input no:''')
                    if j=='no':
                        continue
                    elif j=='yes':
                        break
            g_acc=G*Mass/(R_of_planet**2)
            while True:
                omega_input=int(input('enter tha value of angular velocity of the planet'))
                exp1=int(input('enter the exponent of 10 used in the   :'))
                power1=10**(exp1)
                omega=omega_input*power1
                print('the required angular velocity is:',omega)
                j=input('''If the preview is the desired radius then input yes otherwise input no:''')
                if j=='no':
                    continue
                elif j=='yes':
                    break
            for Q in range(0,91):
                g_actual=g_rotation(math.radians(Q),omega,g_acc,R_of_planet)
                glist.append(g_actual)
            print('g at pole is',g_rotation(math.radians(0),omega,g_acc,R_of_planet))
            gpole=g_rotation(math.radians(0),omega,g_acc,R_of_planet)
            print('g at equater is',g_rotation(math.radians(90),omega,g_acc,R_of_planet))
            geq=g_rotation(math.radians(90),omega,g_acc,R_of_planet)
            r.append(Mass)
            r.append(gpole)
            r.append(geq)
            record.append(r)
            print(str(r))
            file.write('\n'+str(r))
            print("record updated")
            file.close()
            u=input('Do you want to see the graph(yes/no):')
            if u=='yes':
                #plotting the graph using matplotlib module
                #we define the axes of the graph
                plt.plot(glist,color='crimson')
                #we name the axes
                plt.xlabel('angle from the pole')
                plt.ylabel('acceleration due to gravity')
                #displays the graphs
                plt.show()
                an=input('Do you want to save the graph(yes/no):')
                if an=='yes':
                    plt.savefig('grotation.png')
                    print(' file saved')
                    print()
                    print('#1 To the menu')
                    print('#2 Terminate the program')
                    j=int(input('Enter your choice(1/2):'))
                    if j==1:
                        continue
                    if j==2:
                        print('Program terminated')
                        break
                if an=='no':
                    print()
                    print('#1 To the menu')
                    print('#2 Terminate the program')
                    j=int(input('Enter your choice(1/2):'))
                    if j==1:
                        continue
                    if j==2:
                        print('Program terminated')
                        break
                    
            elif u=='no':
                print()
                print('#1 To the menu')
                print('#2 Terminate the program')
                j=int(input('Enter your choice(1/2):'))
                if j==1:
                    continue
                if j==2:
                    print('Program terminated')
                    break
            continue
        elif ch==5:
            print('\t\t ------------------------------------------')
            print('\t\t Gravitational potential energy of a planet')
            print('\t\t ------------------------------------------')
            f=open("potential1.text",'a')
            record=[]
            r=[]
            import matplotlib.pyplot as plt
            graph_list=[]
            print('''
    It is the scalar quantity characteristic of a point in a gravitational field whose
    gradient equals the intensity of the field and equal to the work required to
    move a body of unit mass from given point to a point infinitely remote''')
            print('''
    Consider a thin uniform solid sphere of the radius (R) and mass (M) situated in space.

    Case 1:

    If point ‘P’ lies Inside the uniform solid sphere (r < R): 
    Inside the uniform solid sphere, E = -GMr/R3.
    Using the relation V=-∫ E.dr over a limit of (0 to r).
    The value of gravitational potential is given by,
                     "'V = -GM [(3R2 – r2)/2R2]'"
                     
    Case 2:
    If point ‘P’ lies On the surface of the uniform solid sphere ( r = R ):
    On the surface of a uniform solid sphere, E = -GM/R2. Using the relationV=-∫ E.dr over a limit of (0 to r).
                     "'V = -GM/R.'"

    Case 3:
    If point ‘P’ lies Outside the uniform solid sphere ( r> R): 
    Using the relation over a limit of (0 to r) we get,
                     "'V = -GM/R.'"

    Case 4: Gravitational potential at the centre of the solid sphere is given by,
                     "'V = -3/2 × (GM/R).'"''')
            G=6.67*(10**-11)
            #function to find mass of planet using density
            def mass(radius,rho):
                m_o_planet=rho*(4*3.14*(radius**3))/3
                return m_o_planet
            #we input the values
            while True:
                radius_input=int(input('enter the the significant figures in radius of the planet in m  :'))
                exp1=int(input('enter the exponent of 10 used in the radius  :'))
                power1=10**(exp1)
                R_of_planet=radius_input*power1
                print('the required radius is:',R_of_planet)
                j=input('''If the preview is the desired radius then input yes otherwise input no:''')
                if j=='no':
                    continue
                elif j=='yes':
                    break
            print('if you know the mass of the planet , enter M')
            print('if you know the density of the planet enter D')
            option=input('enter your choice  :')
            if option=='M':
                        #we input significant of mass and the power separately to reduce the work of the user
                while True:
                    mass=int(input('enter the significant figures of mass of the planet  :'))
                    exp=int(input('enter the exponent of 10 used in the mass  :'))
                    power=10**(exp)
                    Mass=mass*power
                    print('the required mass is:',Mass)
                    j=input('''If the preview is the desired mass then input yes otherwise input no:''')
                    if j=='no':
                        continue
                    elif j=='yes':
                        break
            elif option=='D':
                while True:
                    print('''We find the mass of the planet using the Density of the planet using the formula
                        (4*3.14*(Radius**3)/3)*rho''')
                    rho_planet=int(input('enter the density of the planet in kg/m3'))
                    Mass=mass(R_of_planet,rho_planet)
                    print('the required mass is:',Mass)
                    j=input('''If the preview is the desired mass then input yes otherwise input no:''')
                    if j=='no':
                        continue
                    elif j=='yes':
                        break
            while True:
                radius=int(input('enter the significant digits of radius till which field is to calculated'))
                #we input significant of radius and the power separately to reduce the work of the user
                exp1=int(input('enter the exponent of 10 used in the radius  :'))
                power1=10**(exp1)
                r_potential=radius*power1
                print('the  radius is:',r_potential)
                j=input('''If the preview is the desired radius then input yes otherwise input no:''')
                if j=='no':
                    continue
                elif j=='yes':
                    break
            #function to find gravitational potential energy
            def GPE(r,m,R):
                if r>=R:
                    Gravi_potential=(-G*m)/r
                else:
                    Gravi_potential=(-(G*m*(3*(R**2)-(r**2))))/2*(R**2)
                return Gravi_potential
        #finding values of GPE at diff radii
            for i in range(1,r_potential,100000):
                graph_list.append(GPE(i,Mass,R_of_planet))

            print('GPE at surface is',GPE(R_of_planet,Mass,R_of_planet))
            Gsur=GPE(R_of_planet,Mass,R_of_planet)
            r.append(R_of_planet)
            r.append(Mass)
            r.append(Gsur)
            record.append(r)
            print(str(r))
            f.write('\n'+str(r))
            print("record updated")
            f.close()
            
            u=input('Do you want to see the graph(yes/no):')
            if u=='yes':
                #plotting the graph using matplotlib module
            #we define the axes of the graph
                plt.plot(graph_list,color='crimson')
            #we name the axes
                plt.xlabel('radius from center')
                plt.ylabel('Gravitional potential energy')
            #displays the graphs
                plt.show()
                an=input('Do you want to save the graph(yes/no):')
                if an=='yes':
                    plt.savefig('GPE.png')
                    print(' file saved')
                    print()
                    print('#1 To the menu')
                    print('#2 Terminate the program')
                    j=int(input('Enter your choice(1/2):'))
                    if j==1:
                        continue
                    if j==2:
                        print('Program terminated')
                        break
                if an=='no':
                    print()
                    print('#1 To the menu')
                    print('#2 Terminate the program')
                    j=int(input('Enter your choice(1/2):'))
                    if j==1:
                        continue
                    if j==2:
                        print('Program terminated')
                        break
                    
            elif u=='no':
                print()
                print('#1 To the menu')
                print('#2 Terminate the program')
                j=int(input('Enter your choice(1/2):'))
                if j==1:
                    continue
                if j==2:
                    print('Program terminated')
                    break
            continue
        elif ch==6:
            print('\t\t\t -----------------------')
            print('\t\t\t Max height of an object')
            print('\t\t\t -----------------------')

            #max height of a object
            import matplotlib.pyplot as plt
            print('''
    This code uses matplotlib module to make the bar graph of the max heights of a structure that can be formmed by a particular material It uses predefined gravitational
    field of a planet of our solar system along with the elasticity of a particular material to find a maximum height of a structure The stressdue to all the material on
    the top should be less than the critical shearing stress at which the material flow. At the bottom of the object height h ,the force perunit area due to the weight
    of the object is hog where p is the density of the material of the object and g is the acceleration due to the gravity . The materialat the bottom experiences this
    force in the vertical direction and the sides of the  objects are free. Therefore this is not a case of pressure of bulk compression.

    There is the shear component. the shear component streches the object to its sides and thus limits the height of the object

    We can also add our own custom planets and compare it to the planets of our solar system''')
            print()
            while True:
                density=int(input('enter the density of the material  :'))
                print('the density is',density)
                j=input('''If the preview is the desired density then input yes otherwise input no:''')
                if j=='no':
                    continue
                elif j=='yes':
                    break
        #we input significant of elastic limit and the power separately to reduce the work of the user
            while True:
                elastic_=int(input('enter the  the significant figures in elastic limit of the materia  :'))
                exp_=int(input('enter the exponent of 10 used in the radius  :'))
                power_=10**(exp_)
                elastic=elastic_*power_
                print('elastic limit is',elastic)
                j=input('''If the preview is the desired elastic limit then input yes otherwise input no:''')
                if j=='no':
                    continue
                elif j=='yes':
                    break
            print()
            g=[9.8,3.6,8.9,3.8,26.0,11.1,10.7,14.1,275]
            name=['earth','Mercury','Venus','Mars','Jupiter','Saturn','Uranus','Neptune','Sun']
            G=6.67*(10**-11)
            def g_surface(m,r):
                gr=G*m/(r**2)
                print('the acceleration due to gravity at the surface of the planet is ',gr)
                print(''' We find the acceleration due to gravity over the gravitational field of a planet by the help of the formula
                     g=G*mass/R**2''')
                return gr  
            while True:
                ans=input('do u want to input custom planet(y/n)  :')
                if ans=='y':
                    name_planet=input('enter the name of the planet  :')
                    name.append(name_planet)
                    while True:
                        radius=int(input('enter the the significant figures in radius of the planet in m  :'))
                    #we input significant of radius and the power separately to reduce the work of the user
                        exp1=int(input('enter the exponent of 10 used in the radius  :'))
                        power1=10**(exp1)
                        Radius=radius*power1
                        print('the  radius is:',Radius)
                        j=input('''If the preview is the desired radius then input yes otherwise input no:''')
                        if j=='no':
                            continue
                        elif j=='yes':
                            break
                    print('if you know the mass of the planet , enter M')
                    print('if you know the density of the planet enter D')
                    option=input('enter your choice  :')
                    if option=='M':
                    #we input significant of mass and the power separately to reduce the work of the user
                        while True:
                            mass=int(input('enter the significant figures of mass of the planet  :'))
                            exp=int(input('enter the exponent of 10 used in the mass  :'))
                            power=10**(exp)
                            Mass=mass*power
                            print('the required mass is:',Mass)
                            j=input('''If the preview is the desired mass then input yes otherwise input no:''')
                            if j=='no':
                                continue
                            elif j=='yes':
                                break                
                        g_planet=g_surface(Mass,Radius)
                        g.append(g_planet)
                    elif option=='D':
                        print('''We find the mass of the planet using the Density of the planet using the formula
                        (4*3.14*(Radius**3)/3)*rho''')
                        while True:
                            rho=int(input('enter the density of the planet in kg/m3'))
                            Mass=(4*3.14*(Radius**3)/3)*rho
                            print('the required mass is:',Mass)
                            j=input('''If the preview is the desired mass then input yes otherwise input no:''')
                            if j=='no':
                                continue
                            elif j=='yes':
                                break
                        g_planet=g_surface(Mass,Radius)
                        g.append(g_planet)
                    else:
                        print('wrong option')
                        continue
                else:
                    break
            height=[]
            print('''
     We use this formula

     rho*g*h=elastic limit
     h=elasticlimit/(rho*g)

     to find the maximum height of the structures using its elastic city and density''')
            print(density)
            print(g)
            for i in g:

                h=elastic/(density*i)
                height.append(h)
            print()
            u=input('Do you want to see the graph(yes/no):')
            if u=='yes':
                fig = plt.figure()
            #we define the axes of the graph
                ax = fig.add_axes([0,0,1,1])
                ax.bar(name,height)
            #we name the graph
                plt.title('the max height')
            #we name the axes
                plt.xlabel('plantes')
                plt.ylabel('height')
                plt.show() #displays the graphs
                an=input('Do you want to save the graph(yes/no):')
                if an=='yes':
                    plt.savefig('MaxH.png')
                    print(' file saved')
                    print()
                    print('#1 To the menu')
                    print('#2 Terminate the program')
                    j=int(input('Enter your choice(1/2):'))
                    if j==1:
                        continue
                    if j==2:
                        print('Program terminated')
                        break
                if an=='no':
                    print()
                    print('#1 To the menu')
                    print('#2 Terminate the program')
                    j=int(input('Enter your choice(1/2):'))
                    if j==1:
                        continue
                    if j==2:
                        print('Program terminated')
                        break
                    
            elif u=='no':
                print()
                print('#1 To the menu')
                print('#2 Terminate the program')
                j=int(input('Enter your choice(1/2):'))
                if j==1:
                    continue
                if j==2:
                    print('Program terminated')
                    break
            continue
        elif ch==7:
            print('\t\t\t\t ----------------')
            print('\t\t\t\t Orbital Velocity')
            print('\t\t\t\t ----------------')
            import matplotlib.pyplot as plt
            print('''
    Consider a satellite with mass Msat orbiting a central body with a mass of mass MCentral.
    The central body could be a planet, the sun or some other large mass capable of causing sufficient acceleration on a less massive nearby object.
    If the satellite moves in circular motion, then the net centripetal force acting upon this orbiting satellite is given by the relationship

            Fnet = ( Msat • v2 ) / R

    This net centripetal force is the result of the gravitational force that attracts the satellite towards the central body and can be represented as

            Fgrav = ( G • Msat • MCentral ) / R2

    Since Fgrav = Fnet, the above expressions for centripetal force and gravitational force can be set equal to each other. Thus,

            (Msat • v2) / R = (G • Msat • MCentral ) / R2

    Observe that the mass of the satellite is present on both sides of the equation; thus it can be canceled by dividing through by Msat.
    Then both sides of the equation can be multiplied by R, leaving the following equation.

            v2 = (G • MCentral ) / R

    Taking the square root of each side, leaves the following equation for the velocity of a satellite moving about a central body in circular motion

            v=((G • MCentral ) / R)**1/2

    where
    G is 6.673 x 10-11 N•m2/kg2,
    Mcentral is the mass of the central body about which the satellite orbits, and
    R is the radius of orbit for the satellite.
            ''')
            print()
            v_orbital=[]
            G=6.67*(10**-11)
            #we input significant of radius and the power separately to reduce the work of the user
            while True:
                radius_input=int(input('enter the the significant figures in radius of the planet in m  :'))
                exp1=int(input('enter the exponent of 10 used in the radius  :'))
                power1=10**(exp1)
                R_planet=radius_input*power1
                print('the required radius is:',R_planet)
                j=input('''If the preview is the desired radius then input yes otherwise input no:''')
                if j=='no':
                    continue
                elif j=='yes':
                      break
            print()
            #function to find mass of planet using density
            def M(radius,rho):
                m_o_planet=rho*(4*3.14*(radius**3))/3
                return m_o_planet

            print('if you know the mass of the planet , enter M')
            print('if you know the density of the planet enter D')
            option=input('enter your choice  :')
            if option=='M':
                        #we input significant of mass and the power separately to reduce the work of the user
                while True:
                    mass=int(input('enter the significant figures of mass of the planet  :'))
                    exp=int(input('enter the exponent of 10 used in the mass  :'))
                    power=10**(exp)
                    m_planet=mass*power
                    print('the required mass is:',m_planet)
                    j=input('''If the preview is the desired mass then input yes otherwise input no:''')
                    if j=='no':
                         continue
                    elif j=='yes':
                           break                
            elif option=='D':
                while True:
                    print('''We find the mass of the planet using the Density of the planet using the formula
                            (4*3.14*(Radius**3)/3)*rho''')
                    rho_planet=int(input('enter the density of the planet in kg/m3'))
                    m_planet=M(R_planet,rho_planet)
                    print('the required mass is:',m_planet)
                    j=input('''If the preview is the desired mass then input yes otherwise input no:''')
                    if j=='no':
                        continue
                    elif j=='yes':
                        break                
            while True:
                radiusinput=int(input('enter the the significant figures in radius till which the orbit is to be calculated in m  :'))
                exp2=int(input('enter the exponent of 10 used in the radius  :'))
                power2=10**(exp2)
                r_orbit=radiusinput*power2
                print('the required radius is:',r_orbit)
                j=input('''If the preview is the desired radius then input yes otherwise input no:''')
                if j=='no':
                    continue
                elif j=='yes':
                    break                
            def orbit(m,r):
                v=(G*m/r)**0.5
                return v
            for i in range(R_planet,r_orbit,10000):
                vorb=orbit(m_planet,i)
                v_orbital.append(vorb)
            for a in v_orbital:
                 print(a)
            u=input('Do you want to see the graph(yes/no):')
            if u=='yes':
                #plotting the graph using matplotlib module
                #we define the axes of the graph
                plt.plot(v_orbital,color='crimson')
                #we name the axes
                plt.xlabel('radius from center')
                plt.ylabel('orbital velocity')
                #displays the graphs
                plt.show()
                an=input('Do you want to save the graph(yes/no):')
                if an=='yes':
                    plt.savefig('Orbitalvelocity.png')
                    print(' file saved')
                    print()
                    print('#1 To the menu')
                    print('#2 Terminate the program')
                    j=int(input('Enter your choice(1/2):'))
                    if j==1:
                        continue
                    if j==2:
                        print('Program terminated')
                        break
                if an=='no':
                    print()
                    print('#1 To the menu')
                    print('#2 Terminate the program')
                    j=int(input('Enter your choice(1/2):'))
                    if j==1:
                        continue
                    if j==2:
                        print('Program terminated')
                        break
                    
            elif u=='no':
                print()
                print('#1 To the menu')
                print('#2 Terminate the program')
                j=int(input('Enter your choice(1/2):'))
                if j==1:
                    continue
                if j==2:
                    print('Program terminated')
                    break
            continue
        elif ch==8:
            print('\t\t\t\t ----------------')
            print('\t\t\t\t Escape velocity')
            print('\t\t\t\t ----------------')
            #the escape velocity of a planet
            import matplotlib.pyplot as plt
            print('''
    Escape Velocity is referred to as the minimum velocity needed by anybody or object to be projected to overcome the gravitational pull of the planet .
    In other words, the minimum velocity that one requires to escape the gravitational field is escape velocity.

    The formula for escape velocity comprises of a constant, G, which we refer to as the universal gravitational constant.
    The value of G is = 6.673 × 10-11 N . m2 / kg2.
    The unit for escape velocity is meters per second (m/s).

    as energy of the object remains conserved
            -GMm/R+m*v**2/2=0

    thus:
                ve=(2GM/R)**1/2

    Over here:
    ve ---> escape velocity (m/s)
    G  ---> universal gravitational constant (6.673 × 10-11 N . m2 / kg2)
    M  ---> mass of the planet or moon (kg)
    R  ---> radius of the planet or moon (m)''')
            name=['Mercury','Venus','Earth','Moon','Mars','Jupiter','Saturn','Uranus','Neptune',]
            density=[5427,5243,5514,3340,3933,1326,687,1271,1638]
            mass=[0.330,4.87,5.97,0.073,0.642,1898,568,86.8,102]
            diameter=[4879,12104,12756,3475,6792,142984,120536,51118,49528]
            radius=[]
            V_escape=[]
            Actualmass=[]
            for i in diameter:
                r=i*1000*0.5
                radius.append(r)
            for i in mass:
                m=10**24
                M=i*m
                Actualmass.append(M)
            G=6.67*(10**-11)
            while True:
                ans=input('do u want to input custom planet(y/n)  :')
                if ans=='y':
                    name_planet=input('enter the name of the planet  :')
                    name.append(name_planet)
                    while True:
                        radius_input=int(input('enter the the significant figures in radius of the planet in m  :'))
                        #we input significant of radius and the power separately to reduce the work of the user
                        exp1=int(input('enter the exponent of 10 used in the radius  :'))
                        power1=10**(exp1)
                        Radius=radius_input*power1
                        radius.append(Radius)
                        print('the  radius is:',Radius)
                        j=input('''If the preview is the desired radius then input yes otherwise input no:''')
                        if j=='no':
                            continue
                        elif j=='yes':
                            break                
                    print('if you know the mass of the planet , enter M')
                    print('if you know the density of the planet enter D')
                    option=input('enter your choice  :')
                    if option=='M':
                        while True:
                            #we input significant of mass and the power separately to reduce the work of the user
                            mass=int(input('enter the significant figures of mass of the planet  :'))
                            exp=int(input('enter the exponent of 10 used in the mass  :'))
                            power=10**(exp)
                            Mass=mass*power
                            print('the required mass is:',Mass)
                            j=input('''If the preview is the desired mass then input yes otherwise input no:''')
                            if j=='no':
                                continue
                            elif j=='yes':
                                break                
                        Actualmass.append(Mass)
                    elif option=='D':
                        print('''We find the mass of the planet using the Density of the planet using the formula
                        (4*3.14*(Radius**3)/3)*rho''')
                        while True:
                            rho=int(input('enter the density of the planet in kg/m3'))
                            Mass=(4*3.14*(Radius**3)/3)*rho
                            print('the required mass is:',Mass)
                            j=input('''If the preview is the desired mass then input yes otherwise input no:''')
                            if j=='no':
                                continue
                            elif j=='yes':
                                break                
                        Actualmass.append(Mass)
                    else:
                        print('wrong option')
                        continue
                else:
                    break
            print(Actualmass)
            print(radius)
            print(name)
            def escape_v(M,R):
                ve=(2*G*M/R)**(1/2)
                return ve
            length=len(name)
            for i in range (0,length):
                v=escape_v(Actualmass[i],radius[i])
                V_escape.append(v)
                print(V_escape[i],'for',name[i])
            u=input('Do you want to see the graph(yes/no):')
            if u=='yes':
                fig = plt.figure()
            #we define the axes of the graph
                ax = fig.add_axes([0,0,1,1])
                ax.bar(name,V_escape)
            #we name the graph
                plt.title('the escape v')
                #we name the axes
                plt.xlabel('plantes')
                plt.ylabel('v')
                plt.show() #displays the graphs
                an=input('Do you want to save the graph(yes/no):')
                if an=='yes':
                    plt.savefig('grotation.png')
                    print(' file saved')
                    print()
                    print('#1 To the menu')
                    print('#2 Terminate the program')
                    j=int(input('Enter your choice(1/2):'))
                    if j==1:
                        continue
                    if j==2:
                        print('Program terminated')
                        break
                if an=='no':
                    print()
                    print('#1 To the menu')
                    print('#2 Terminate the program')
                    j=int(input('Enter your choice(1/2):'))
                    if j==1:
                        continue
                    if j==2:
                        print('Program terminated')
                        break
                    
            elif u=='no':
                print()
                print('#1 To the menu')
                print('#2 Terminate the program')
                j=int(input('Enter your choice(1/2):'))
                if j==1:
                    continue
                if j==2:
                    print('Program terminated')
                    break
            continue
        elif ch==9:
            print('\t\t\t\t -------------------------')
            print('\t\t\t\t Time period of a satalite')
            print('\t\t\t\t -------------------------')
            f=open("gravity1.text",'a')
            record=[]
            r=[]

            print('''
    The period of a satellite (T) and the mean distance from the central body (R) are related by the following equation:

    T^2/R^3=4*pi^2/G*M

    where
    T is the period of the satellite,
    R is the average radius of orbit for the satellite (distance from center of central planet), and
    G is 6.673 x 10-11 N•m2/kg2.

    There is an important concept evident in all three of these equations -
    the period, speed and the acceleration of an orbiting satellite are not dependent upon the mass of the satellite.

    None of these three equations has the variable Msatellite in them.
    The period, speed and acceleration of a satellite are only dependent upon the radius of orbit and the mass of the central body that the satellite is orbiting.
    Just as in the case of the motion of projectiles on earth, the mass of the projectile has no effect upon the acceleration towards the earth and the speed at
    any instant. When air resistance is negligible and only gravity is present, the mass of the moving object becomes a non-factor.
    Such is the case of orbiting satellites.''')
            print()
            G=6.67*(10**-11)
            #we input significant of radius and the power separately to reduce the work of the user
            radius_input=int(input('enter the the significant figures in radius of the planet in m  :'))
            exp1=int(input('enter the exponent of 10 used in the radius  :'))
            power1=10**(exp1)
            R_planet=radius_input*power1
            r.append(R_planet)
            print('the required radius is:',R_planet)
            print()
            #function to find mass of planet using density
            def M(radius,rho):
                m_o_planet=rho*(4*3.14*(radius**3))/3
                return m_o_planet       
     
            print('if you know the mass of the planet , enter M')
            print('if you know the density of the planet enter p')
            option=input('enter your choice  :')
            if option=='M':
                        #we input significant of mass and the power separately to reduce the work of the user
                mass=int(input('enter the significant figures of mass of the planet  :'))
                exp=int(input('enter the exponent of 10 used in the mass  :'))
                power=10**(exp)
                m_planet=mass*power
                print('the required mass is:',m_planet)
            elif option=='p':
                print('''We find the mass of the planet using the Density of the planet using the formula
    (4*3.14*(Radius**3)/3)*rho''')
                rho_planet=int(input('enter the density of the planet in kg/m3'))
                m_planet=M(R_planet,rho_planet)
                r.append(m_planet)
                print('the required mass is:',m_planet)
            radiusinput=int(input('enter the the significant figures in radius till which the orbit is to be calculated in m  :'))
            exp2=int(input('enter the exponent of 10 used in the radius  :'))
            power2=10**(exp2)
            r_orbit=radiusinput*power2
            print('the required radius is:',r_orbit)
            def time(m,r):
                t=(((r**3)*4*3.14*3.14)/(G*m))**0.5
                return t
            r.append(time (m_planet,r_orbit))
            record.append(r)
            f.write('\n'+str(r))
            print('time period is ',time(m_planet,r_orbit))
            print("record updated")
            f.close()
            
        elif ch==10:
            import math
            print()
            print('\t\t\t -----------------------')
            print('\t\t\t Altitude of Kármán Line')
            print('\t\t\t -----------------------')
            print('''The Kármán line is the altitude where space begins. It is 100 km(62 miles) high.
    It represents the boundary between Earth's atmosphere and outer space.
    The Kármán line is named after Theodore von Kármán, a Hungarian-American engineer and physicist.''')
            print()
            print('''At this line atmosphere becomes too thin to support aeronautical flight. An aircraft
    at this altitude would have to travel faster than orbital velocity to obtain enough
    lift to support itself.''')
            print()
            print('''There is a sudden increase in temperature of the atmosphere and solar radiation just
    below the line. This places the line within the greater thermosphere.''')
            print()
            print('''Kármán line is usually calculated for Earth but it is also calculated for other
    planets during interplanetary missions.''')
            print()
            k=1.38*(10**-23)
            G=6.67*(10**-11)
            while True:
                r=float(input('Enter the radius of the planet:'))
                print('Preview of the radius:')
                print(r,'m',sep='')
                n=input('If the preview is the desired radius then input yes otherwise input no:')
                if n=='yes':
                    R=r
                    break
                elif n=='no':
                    continue
            print()
            t=input('''If the mass of the planet is known then input M
    If the average density of the planet is known then input p:''')
            if t=='M':
                print()
                print('''# Please enter the mass up to second decimal place''')
                while True:
                    print()
                    M=float(input('Enter the mass:'))
                    e=int(input('Enter the power of 10:'))
                    print('Preview of the mass:')
                    z=M*(10**e)
                    print(z,'kg',sep='')
                    ch=input('If the preview is the desired mass then input yes otherwise input no:')
                    if ch=='yes':
                        z=M*(10**e)
                        break
                    elif ch=='no':
                        continue
                g=(G*z)/R**2
            elif t=='p':
                p=int(input('Enter the density(kg/m^3):'))
                z=(4/3)*3.14*(R**3)*p
            g=(G*z)/R**2
            an=6.023*(10**26)
            esp=14.7
            ep=32000
            print()
            P=float(input('Enter the atmospheric pressure on the planet(psi):'))
            print()
            T=int(input('Enter the surface temperature of the planet(K):'))
            print()
            percentage=int(input('Enter the percentage of the major constituent(compound) of the atmosphere(%):'))
            per=100-percentage
            m=int(input('Enter the molar mass of the major compound(amu):'))
            mo=m/an
            A=(k*T)/(mo*g)
            x=(ep*P)/esp
            y=(x*m)/28
            h=A*(math.log(y))
            print()
            j=input('''If the molar mass of other atmospheric constituents is less then
    the major, then input yes otherwise input no:''')
            if j=='yes':
                hf=h+((per/100)*h)
            elif j=='no':
                hf=h-((per/100)*h)
            print()
            print('The altitude of the Kármán line for the given details is',round((hf/1000),2),'km.' )
            if per>=20:
                print()
                print('''# There maybe be some error in the value of altitude(+/- 30% approx)
    because of higher percentage of other gases and change in temperature with altitude.''')
            elif per<20:
                print()
                print('''# There maybe be some error in the value of altitude(+/- 10% approx)
    because the change in temperature with altitude is ignored.''')
            print()
            print('#1 To the menu')
            print('#2 Terminate the program')
            j=int(input('Enter your choice(1/2):'))
            if j==1:
                continue
            if j==2:
                print('Program terminated')
                break
            continue
        elif ch==11:
            print()
            print('\t\t\t\t ----------------------')
            print('\t\t\t\t Gravitational Redshift')
            print('\t\t\t\t ----------------------')
            print()
            print('''
    Einstein’s theory of general relativity predicts that the wavelength of electromagnetic
    radiation will lengthen as it climbs out of a gravitational well.''')
            print()
            print('''
    Photons must expend energy to escape, but at the same time must always travel at the speed of light,
    so this energy must be lost through a change of frequency rather than a change in speed.
    If the energy of the photon decreases, the frequency also decreases.

    This corresponds to an increase in the wavelength of the photon, or a shift to the red end of the
    electromagnetic spectrum. Thus the name Gravitational Redshift.''')
            print()
            print()
            G=6.67*(10**-11)
            while True:
                r=float(input('Enter the radius of the planet:'))
                print('Preview of the radius:')
                print(r,'m',sep='')
                n=input('If the preview is the desired radius then input yes otherwise input no:')
                print()
                if n=='yes':
                    R=r
                    break
                elif n=='no':
                    continue
            t=input('''
    If the mass of the planet is known then input M
    If the average density of the planet is known then input p:''')
            if t=='M':
                print()
                print('''# Please enter the mass up to second decimal place''')
                while True:
                    print()
                    M=float(input('Enter the mass:'))
                    e=int(input('Enter the power of 10:'))
                    print('Preview of the mass:')
                    z=M*(10**e)
                    print(z,'kg',sep='')
                    ch=input('If the preview is the desired mass then input yes otherwise input no:')
                    if ch=='yes':
                        z=M*(10**20)
                        break
                    elif ch=='no':
                        continue
                g=(G*z)/R**2
            elif t=='p':
                print()
                p=int(input('Enter the density(kg/m^3):'))
                z=(4/3)*3.14*(R**3)*p
            g=(G*z)/R**2
            print()
            print('''The frequency of the radiated electromagnetic was decreases with an
    increase in its distance fom the source of the gravitational field.''')
            print()
            h=int(input('Enter the height from the surface of the planet(in metres):'))
            print()
            print('For reference-')
            print('Hz(10^0)')
            print('kHz(10^3)')
            print('MHz(10^6)')
            print('GHz(10^9)')
            print('THz(10^12)')
            print('PHz(10^15)')
            print('EHz(10^16)')
            print()
            l=[10,10**3,10**6,10**9,10**12]
            f=float(input('Enter the frequency of the wave upto 2 decimal places:'))
            y=int(input('Enter the power of 10:'))
            x=10**y
            fi=f*x
            c=299792458
            x=[]
            for a in range(0,h+1,(h//50)):
                ff=((g*a*fi)/(c**2))+fi
                print(ff,'Hz')
            print()
            print('Initial frequency=',fi,'Hz')
            print('Final frequency=',ff,'Hz')
            print()
            print('''The concept of Gravitational Redshift is very useful in determining
    the distance and velocity of far away celestial bodies from earth.''')
            print()
            print('#1 To the menu')
            print('#2 Terminate the program')
            j=int(input('Enter your choice(1/2):'))
            if j==1:
                continue
            if j==2:
                print()
                print('Program terminated')
                break

        elif ch==12:
            print()
            print('''
    \t\t\t ------------------
    \t\t\t Program Terminated
    \t\t\t ------------------''')
            print()
            print('''
    \t\t\t  ----------------
    \t\t\t  Thanks for using
    \t\t\t  ----------------
    ''')
            break
except ImportError:
    print('download kar na')
