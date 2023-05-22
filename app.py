import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

g = 9.81  # Accélération due à la gravité


def compute_ball_coordinates(v0, angle, t):
    radian_angle = np.deg2rad(angle)
    x = v0 * np.cos(radian_angle) * t
    y = v0 * np.sin(radian_angle) * t - 0.5 * g * t**2
    return x, y


st.title("Courbe de ballon de rugby")
v0 = st.number_input("Vitesse initiale (en m/s)", min_value=0.0, step=1.0, value=10.0)
angle = st.number_input("Angle initial (en degrés)", min_value=0.0, max_value=90.0, step=1.0, value=45.0)

t = np.linspace(0, 2 * v0 * np.sin(np.deg2rad(angle)) / g, num=100)  # Intervalles de temps
x, y = compute_ball_coordinates(v0, angle, t)
plt.plot(x, y)
plt.xlabel("Distance (m)")
plt.ylabel("Hauteur (m)")
plt.title("Trajectoire du ballon de rugby")
st.pyplot()