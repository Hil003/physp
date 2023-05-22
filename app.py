import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# Chemin d'accès à l'image
image_path = "Desktop/physp/"

# Code HTML/CSS pour l'image de fond
html_code = f'''
<style>
body {{
    background-image: url("{image_path}");
    background-size: cover;
}}
</style>
'''

# Affichage de l'image de fond
st.markdown(html_code, unsafe_allow_html=True)

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



"""
sur streamlit
"""
st.title("Courbe de ballon de rugby")
v0 = st.number_input("Vitesse initiale (en m/s)", min_value=0.0, step=1.0, value=10.0)
angle = st.number_input("Angle initial (en degrés)", min_value=0.0, max_value=90.0, step=1.0, value=45.0)

tflight, xmax, ymax = search_t_flight(v0, angle)

print(type(float(tflight)))
tfinal = st.slider("Temps (en secondes)", min_value=0.0, max_value=float(tflight), step=0.01, value=float(tflight/2))



y = 0

"""
debugage

v0 = 10
angle = 45
"""


print(tflight, xmax)
X, Y = [], []
t = 0
step = 0.001

while t < tfinal:
    t, x, y = compute_ball_coordinates(v0, angle, t)
    X.append(x)
    Y.append(y)
    t += step

fig, ax = plt.subplots()
plt.plot(X, Y)
plt.ylim(-ymax*0.05, ymax*1.1)
plt.xlim(-xmax*0.05, xmax*1.1)
plt.xlabel("Distance (m)")
plt.ylabel("Hauteur (m)")
plt.title("Trajectoire du ballon de rugby en fonction du temps")
# plt.show()
st.pyplot(fig)
