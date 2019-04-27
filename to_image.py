from PIL import Image
import numpy
import csv


def main():
    fname = "image-0.dat"
    img = []

    with open(fname) as f:
        for line in f:
            d = line.split(" ")
            e = [float(i) for i in d]
            img.append(e)
    arr_img = numpy.array(img, ndmin=2) * 255
    img_pil = Image.fromarray(arr_img)
    img_pil.show()


main()