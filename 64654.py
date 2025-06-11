import math
import numpy as np
import matplotlib.pyplot as plt

R = 4  # Crystal cell length set to 4, can be divided into 4 parts

class K3C60:
    def __init__(self):
        self.m = 0
        self.l = 0
        self.i = 0
        self.j = 0
        self.k = 0
        self.flag = 0
        self.n1 = 0
        self.n2 = 0
        self.n3 = 0
        self.x = 0
        self.r = 0.0
        self.a = 0.0
        self.ratio = 1.0
        self.aver = 0.0
        self.a_array = []
        self.cnt = 0
        self.avg_array = []
        self.avg_cnt = []

    def convertflag_k(self, x):
        return 1 if x % 2 == 1 else -3

    def convertflag_c60(self, x):
        return 9 if x % 2 == 1 else -3

    def result_k1(self, n):
        for l in range(2, (n // 2) + 1):
            self.m = 2 * l
            self.cnt += 1
            if l % 2 == 0:
                # Calculate interaction energy between c60 and k1
                if l % 4 == 0:
                    self.flag = 1
                    self.x = 0
                else:
                    self.flag = -3
                    self.x = 1

                for i in range(-l, l + 1, 2):
                    self.n1 = i
                    for j in range(-l, l + 1, 2):
                        self.n2 = j
                        for k in range(-l, l + 1, 2):
                            self.n3 = k
                            if self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3 == 0:
                                self.flag = self.convertflag_k(self.x)
                                self.x += 1
                                continue
                            self.r = math.sqrt(self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3)
                            if self.n1 == l or self.n1 == -l:
                                self.ratio /= 2
                            if self.n2 == l or self.n2 == -l:
                                self.ratio /= 2
                            if self.n3 == l or self.n3 == -l:
                                self.ratio /= 2
                            self.a += -self.flag * self.ratio * R / self.r
                            self.flag = self.convertflag_k(self.x)
                            self.x += 1
                            self.ratio = 1.0

                # Calculate interaction with k2
                for i in range(-l + 1, l, 2):
                    self.n1 = i
                    for j in range(-l + 1, l, 2):
                        self.n2 = j
                        for k in range(-l + 1, l, 2):
                            self.n3 = k
                            self.r = math.sqrt(self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3)
                            self.a += -1 * R / self.r
            else:
                # Calculate interaction with k2
                for i in range(-l, l + 1, 2):
                    self.n1 = i
                    for j in range(-l, l + 1, 2):
                        self.n2 = j
                        for k in range(-l, l + 1, 2):
                            self.n3 = k
                            self.r = math.sqrt(self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3)
                            if self.n1 == l or self.n1 == -l:
                                self.ratio /= 2
                            if self.n2 == l or self.n2 == -l:
                                self.ratio /= 2
                            if self.n3 == l or self.n3 == -l:
                                self.ratio /= 2
                            self.a += -1 * self.ratio * R / self.r
                            self.ratio = 1.0

                # Calculate interaction with k1, c60
                if (l - 1) % 4 == 0:
                    self.flag = 1
                    self.x = 0
                else:
                    self.flag = -3
                    self.x = 1

                for i in range(-l + 1, l, 2):
                    self.n1 = i
                    for j in range(-l + 1, l, 2):
                        self.n2 = j
                        for k in range(-l + 1, l, 2):
                            self.n3 = k
                            if self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3 == 0:
                                self.flag = self.convertflag_k(self.x)
                                self.x += 1
                                continue
                            self.r = math.sqrt(self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3)
                            self.a += -self.flag * R / self.r
                            self.flag = self.convertflag_k(self.x)
                            self.x += 1
                            self.ratio = 1.0

            self.a_array.append(self.a)
            print(f"m={self.m}\ta={self.a}")
            if self.cnt % 4 == 0:
                self.aver = np.mean(self.a_array)
                self.avg_array.append(self.aver)
                self.avg_cnt.append(self.cnt // 4)
                print(f"K1 result's average = {self.aver}")
            self.a = 0.0

        print()
        print(f"K1 final result's = {self.aver}")
        print()

        # Plotting
        plt.figure(figsize=(8, 6))
        plt.plot(self.avg_cnt, self.avg_array, 'b-o', label='K1 Average',color='black')
        plt.title('K1 Average vs. Average Count')
        plt.xlabel('Average Count')
        plt.ylabel('Average Value')
        plt.grid(True)
        plt.legend()
        plt.savefig('k1_average_plot.png')
        plt.show()
    def result_c60(self, n):
        for l in range(2, (n // 2) + 1):
            self.m = 2 * l
            self.cnt += 1
            if l % 2 == 0:
                # Calculate interaction energy between c60 and k1
                if l % 4 == 0:
                    self.flag = 9
                    self.x = 0
                else:
                    self.flag = -3
                    self.x = 1

                for i in range(-l, l + 1, 2):
                    self.n1 = i
                    for j in range(-l, l + 1, 2):
                        self.n2 = j
                        for k in range(-l, l + 1, 2):
                            self.n3 = k
                            if self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3 == 0:
                                self.flag = self.convertflag_c60(self.x)
                                self.x += 1
                                continue
                            self.r = math.sqrt(self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3)
                            if self.n1 == l or self.n1 == -l:
                                self.ratio /= 2
                            if self.n2 == l or self.n2 == -l:
                                self.ratio /= 2
                            if self.n3 == l or self.n3 == -l:
                                self.ratio /= 2
                            self.a += -self.flag * self.ratio * R / self.r
                            self.flag = self.convertflag_c60(self.x)
                            self.x += 1
                            self.ratio = 1.0

                # Calculate interaction with k2
                for i in range(-l + 1, l, 2):
                    self.n1 = i
                    for j in range(-l + 1, l, 2):
                        self.n2 = j
                        for k in range(-l + 1, l, 2):
                            self.n3 = k
                            self.r = math.sqrt(self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3)
                            self.a += 3 * R / self.r
            else:
                # Calculate interaction with k2
                for i in range(-l, l + 1, 2):
                    self.n1 = i
                    for j in range(-l, l + 1, 2):
                        self.n2 = j
                        for k in range(-l, l + 1, 2):
                            self.n3 = k
                            self.r = math.sqrt(self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3)
                            if self.n1 == l or self.n1 == -l:
                                self.ratio /= 2
                            if self.n2 == l or self.n2 == -l:
                                self.ratio /= 2
                            if self.n3 == l or self.n3 == -l:
                                self.ratio /= 2
                            self.a += 3 * self.ratio * R / self.r
                            self.ratio = 1.0

                # Calculate interaction with k1, c60
                if (l - 1) % 4 == 0:
                    self.flag = 9
                    self.x = 0
                else:
                    self.flag = -3
                    self.x = 1

                for i in range(-l + 1, l, 2):
                    self.n1 = i
                    for j in range(-l + 1, l, 2):
                        self.n2 = j
                        for k in range(-l + 1, l, 2):
                            self.n3 = k
                            if self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3 == 0:
                                self.flag = self.convertflag_c60(self.x)
                                self.x += 1
                                continue
                            self.r = math.sqrt(self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3)
                            self.a += -self.flag * R / self.r
                            self.flag = self.convertflag_c60(self.x)
                            self.x += 1

            self.a_array.append(self.a)
            print(f"m={self.m}\ta={self.a}")
            if self.cnt % 4 == 0:
                self.aver = np.mean(self.a_array)
                self.avg_array.append(self.aver)
                self.avg_cnt.append(self.cnt // 4)
                print(f"C60 result's average = {self.aver}")
                print(f"Debug: a_array={self.a_array}, avg_array={self.avg_array}, avg_cnt={self.avg_cnt}")
            self.a = 0.0  # Reset a for next iteration

        print()
        print(f"C60 final result's = {self.aver}")
        print()

        # Plotting
        plt.figure(figsize=(8, 6))
        if self.avg_cnt and self.avg_array:  # Check if data exists
            plt.plot(self.avg_cnt, self.avg_array, 'r-o', label='C60 Average')
            plt.title('C60 Average vs. Average Count')
            plt.xlabel('Average Count')
            plt.ylabel('Average Value')
            plt.grid(True)
            plt.legend()
            plt.savefig('c60_average_plot.png')
            plt.show()  # Ensure graph is displayed
        else:
            print("Warning: No data to plot for C60 (avg_cnt or avg_array is empty)")

    def result_k2(self, n):
        for l in range(2, (n // 2) + 1):
            self.m = 2 * l
            self.cnt += 1
            if l % 2 == 0:
                for i in range(-l, l + 1, 2):
                    self.n1 = i
                    for j in range(-l, l + 1, 2):
                        self.n2 = j
                        for k in range(-l, l + 1, 2):
                            self.n3 = k
                            if self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3 == 0:
                                continue
                            self.r = math.sqrt(self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3)
                            if self.n1 == l or self.n1 == -l:
                                self.ratio /= 2
                            if self.n2 == l or self.n2 == -l:
                                self.ratio /= 2
                            if self.n3 == l or self.n3 == -l:
                                self.ratio /= 2
                            self.a += -1 * self.ratio * R / self.r
                            self.flag = self.convertflag_k(self.x)
                            self.x += 1
                            self.ratio = 1.0

                # Calculate interaction with k1, c
                for i in range(-l + 1, l, 2):
                    self.n1 = i
                    for j in range(-l + 1, l, 2):
                        self.n2 = j
                        for k in range(-l + 1, l, 2):
                            self.n3 = k
                            self.r = math.sqrt(self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3)
                            self.flag = self.convertflag_k(self.x)
                            self.x += 1
                            self.a += -self.flag * R / self.r
            else:
                # Calculate K1 and C60 interaction
                for i in range(-l, l + 1, 2):
                    self.n1 = i
                    for j in range(-l, l + 1, 2):
                        self.n2 = j
                        for k in range(-l, l + 1, 2):
                            self.n3 = k
                            self.r = math.sqrt(self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3)
                            if self.n1 == l or self.n1 == -l:
                                self.ratio /= 2
                            if self.n2 == l or self.n2 == -l:
                                self.ratio /= 2
                            if self.n3 == l or self.n3 == -l:
                                self.ratio /= 2
                            self.flag = self.convertflag_k(self.x)
                            self.x += 1
                            self.a += -self.flag * self.ratio * R / self.r
                            self.ratio = 1.0

                # Calculate interaction with k2
                for i in range(-l + 1, l, 2):
                    self.n1 = i
                    for j in range(-l + 1, l, 2):
                        self.n2 = j
                        for k in range(-l + 1, l, 2):
                            self.n3 = k
                            if self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3 == 0:
                                continue
                            self.r = math.sqrt(self.n1 * self.n1 + self.n2 * self.n2 + self.n3 * self.n3)
                            self.a += -1 * R / self.r

            self.a_array.append(self.a)
            print(f"m={self.m}\ta={self.a}")
            if self.cnt % 4 == 0:
                self.aver = np.mean(self.a_array)
                self.avg_array.append(self.aver)
                self.avg_cnt.append(self.cnt // 4)
                print(f"K2 result's average = {self.aver}")
            self.a = 0.0

        print()
        print(f"K2 final result's = {self.aver}")
        print()

        # Plotting
        plt.figure(figsize=(8, 6))
        plt.plot(self.avg_cnt, self.avg_array, 'g-o', label='K2 Average')
        plt.title('K2 Average vs. Average Count')
        plt.xlabel('Average Count')
        plt.ylabel('Average Value')
        plt.grid(True)
        plt.legend()
        plt.savefig('k2_average_plot.png')
        plt.show()
def main():
    k1 = K3C60()
    c60 = K3C60()
    k2 = K3C60()

    print("K1:")
    n = int(input("please enter the width of solid(n): "))
    k1.result_k1(n)

    print("C60:")
    n = int(input("please enter the width of solid(n): "))
    c60.result_c60(n)

    print("K2:")
    n = int(input("please enter the width of solid(n): "))
    k2.result_k2(n)

    input("Press Enter to continue...")

if __name__ == "__main__":
    main()