import math



def cnt_radius(m,n):
    c_c = 1.42
    diameter = c_c / math.pi
    diameter = diameter * math.sqrt(3*(m*m+m*n+n*n))

    radius = diameter / 2.0
    print(n,m)
    print("diameter ", str(diameter))
    print("radius " + str(radius))


    return radius

if __name__ == "__main__":
    m = float(input("m:"))
    n = float(input("n:"))


    cnt_radius(m,n)

