import matplotlib.pyplot as plt
from ..constants import get_colors

def add_legend_panels(Obj, i_Emax, D_i_c_):
    fig, ax = plt.subplots()
    ax.set_xlim(0, 200) 
    ax.set_ylim(0, 200)
    colors = get_colors()
    x_a = -10
    dx = 200 / 6
    i =-10 
    for j in range(2,29): 
        if i_Emax[j, 0] == 0: continue
        particle = str(Obj[j][9])
        plt.text(x_a, i, Obj[j][0], fontsize=8)
        plt.text(x_a + 0.7* dx, i, particle, color=colors[j], fontsize=8)
        plt.text(x_a + 1.7 * dx, i, Obj[j][1], fontsize=8)
        plt.text(x_a + 3.2 * dx, i, Obj[j][5], fontsize=8)        
        plt.text(x_a + 4.5 * dx, i, i_Emax[j, 0], fontsize=8)
        plt.text(x_a + 5.5 * dx, i, D_i_c_[j], fontsize=8)
      #  i += 180    # für 333           # 100 für 222
        i += 9
    plt.text(x_a + 1.7 * dx, i, " mass ", fontsize=8)
    plt.text(x_a + 3.2 * dx, i, "average lifespan", fontsize=8)    
    plt.text(x_a + 4.5 * dx, i, "  ∆i  ", fontsize=8)
    plt.text(x_a + 5.5 * dx, i, " ∆i/(2pi) ", fontsize=8)
 
    plt.axis('off')

    return()
