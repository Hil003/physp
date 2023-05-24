import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

g = 9.81  # Accélération due à la gravité

def compute_ball_coordinates(v0, angle, t):
    radian_angle = np.deg2rad(angle)
    x = v0 * np.cos(radian_angle) * t
    y = v0 * np.sin(radian_angle) * t - 0.5 * g * t**2
    return t, x, y


def search_t_flight(v0, angle):
    g = 9.81  # Accélération due à la gravité
    radian_angle = np.radians(angle)
    t_flight = (2 * v0 * np.sin(radian_angle)) / g
    d_max = (v0 ** 2 * np.sin(2 * radian_angle)) / g
    h_max = (v0**2 * np.sin(radian_angle)**2) / (2 * g)

    return t_flight, d_max, h_max





def trace_sans_frott():

    """
    sur streamlit
    """
    st.title("Courbe de ballon de rugby")
    v0 = st.number_input("Vitesse initiale (en m/s)", min_value=0.0, step=1.0, value=10.0)
    angle = st.number_input("Angle initial (en degrés)", min_value=0.0, max_value=90.0, step=1.0, value=45.0)

    tflight, xmax, ymax = search_t_flight(v0, angle)

    print(type(float(tflight)))
    tfinal = st.slider("Temps (en secondes)", min_value=0.0, max_value=float(tflight), step=0.01, value=float(tflight/2))

    """
    debugage

    v0 = 10
    angle = 45
    """
    X, Y = [], []
    t = 0
    step = 0.001

    while t < tfinal:
        t, x, y = compute_ball_coordinates(v0, angle, t)
        X.append(x)
        Y.append(y)
        t += step

    lim = max(xmax, ymax)

    fig, ax = plt.subplots()
    plt.plot(X, Y)
    plt.ylim(-lim*0.05, lim*1.1)
    plt.xlim(-lim*0.05, lim*1.1)
    plt.xlabel("Distance (m)")
    plt.ylabel("Hauteur (m)")
    plt.title("Trajectoire du ballon de rugby en fonction du temps")
    # plt.show()
    st.pyplot(fig)
    return 0