import matplotlib.pyplot as plt

f = open('../SRCone/SRCone/contact_angle.record','r')
x=[]
y=[]

for line in  f.readlines():
    l = line.strip().split(",")
    x.append(float(l[0]))
    y.append(float(l[1]))
#x = [0.0, -0.05, -0.1, -0.15, -0.2, -0.25, -0.3, -0.35, -0.4, -0.45 ,-0.5,-0.6,-0.7,-0.8,-0.9]
#y = [180.0, 180.0, 180.0, 180.0, 160.148, 149.247, 136.918, 129.445, 121.605, 113.45,105.048,92.2466,69.1846,50.9267,35.4893]
f.close()
plt.plot(x,y,marker="o")
#plt.scatter(x,y)
plt.show()