from Exercise2 import plt, run_implicitEuler, imagepath

def f3(y, t):
    """Get the slope of a known segment around point (y_0, t_0)"""
    f = 1/(y^2 +1)
    return f

def df_dy3(y):
    """Gets the derivative of the slope as a function of y(t) at point (y)"""
    df_dy = 2*y*(-1/(y**2 + 1)**2)
    return df_dy


#Call and plot the functions
#Initialise
y_0 = 1
t_0 = 0
err_0 = 0.1
y_1_guess = 1.5 #this is an initial guess

#parameters
steps = 10

#prepare plots
fig2, (ax3, ax4) = plt.subplots(2,1)

#Run for 3 time steps
linecolours = ["orange", "green"]
errcolours = ["orange", "lightgreen"]
for Dt in (2,3):
    lsty_imp, lstf_imp, lstt_imp, lsterr_imp = run_implicitEuler(df_dy3, y_0, t_0, err_0, y_1_guess, Dt, steps)

    # plot y VS t
    ax3.plot(lstt_imp/Dt, lsty_imp, color = linecolours[Dt-2], label = f"time step = {Dt}")
    ax3.scatter(lstt_imp/Dt, lsty_imp, marker ="+", color = "black")

    # Plot the error
    ax3.errorbar(lstt_imp/Dt, lsty_imp, abs(lsterr_imp/max(lsterr_imp)), fmt = "None", color = errcolours[Dt-2], elinewidth = 0.3, capsize = 3)

    # Plot the slope
    ax4.plot(lstt_imp/Dt, lstf_imp, label = f"time step = {Dt}")

#Finish the plots
ax3.set_xlabel("t")
ax3.set_ylabel("y(t)")
ax3.legend()
ax3.grid(alpha=0.5)

ax4.set_xlabel("t")
ax4.set_ylabel("f(y,t)")
ax4.legend()
ax4.grid(alpha=0.5)
fig2.savefig(imagepath+"ex3_implicitEuler.pdf", dpi=300)


# Plot Dt = 1 separately because it obscures the other ones
fig3, ax5 = plt.subplots(1,1, figsize=(10,4))

Dt = 1
lsty_imp, lstf_imp, lstt_imp, lsterr_imp = run_implicitEuler(df_dy3, y_0, t_0, err_0, y_1_guess, Dt, steps)

# plot y VS t
ax5.plot(lstt_imp/Dt, lsty_imp, color = "blue", label = f"time step = {Dt}")
ax5.scatter(lstt_imp/Dt, lsty_imp, marker = "+", color = "black")

# Plot the error
ax5.errorbar(lstt_imp/Dt, lsty_imp, abs(lsterr_imp/max(lsterr_imp)), fmt = "None", color = "lightblue", elinewidth = 0.3, capsize = 3)

#Finish the plots
ax5.set_xlabel("t")
ax5.set_ylabel("y(t)")
ax5.legend()
ax5.grid(alpha=0.5)

fig3.savefig(imagepath+"ex3_implicitEulerDt1.pdf", dpi=300)
