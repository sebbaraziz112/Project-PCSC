
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread

# Charger l’image (grayscale ou couleur)
from PIL import Image

img = Image.open("../ressources/ImageFolder/sine_diag.png")
img.save("../ressources/ImageFolder/bird_fixed.png")
img = imread("../ressources/ImageFolder/bird_fixed.png")

# Si l’image est couleur → convertir en niveaux de gris
if img.ndim == 3:
    img = img.mean(axis=2)

# FFT 2D
F = np.fft.fft2(img)

# Décalage du spectre (centre les basses fréquences)
Fshift = np.fft.fftshift(F)

# Magnitude (amplitude)
magnitude = np.abs(Fshift)

# Log-scale pour affichage
magnitude_log = np.log1p(magnitude)

# Affichage
plt.figure(figsize=(6,6))
plt.imshow(magnitude_log, cmap='gray')
plt.title("Amplitude du spectre (log)")
plt.axis('off')
plt.show()
