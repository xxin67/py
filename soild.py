import math

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
        self.R = 4  # 晶胞长度

    def _convertflag_k(self, x):
        return 1 if x % 2 == 1 else -3

    def _convertflag_c60(self, x):
        return 9 if x % 2 == 1 else -3

    def result_k1(self, n):
        for l in range(2, n // 2 + 1):
            self.l = l
            self.m = 2 * l
            self.cnt += 1
            if l % 2 == 0:
                # 处理C60与K1面
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
                            if self.n1**2 + self.n2**2 + self.n3**2 == 0:
                                self.flag = self._convertflag_k(self.x)
                                self.x += 1
                                continue
                            self.r = math.sqrt(self.n1**2 + self.n2**2 + self.n3**2)
                            self.ratio = 1.0
                            if self.n1 == l or self.n1 == -l:
                                self.ratio /= 2
                            if self.n2 == l or self.n2 == -l:
                                self.ratio /= 2
                            if self.n3 == l or self.n3 == -l:
                                self.ratio /= 2
                            self.a += -self.flag * self.ratio * self.R / self.r
                            self.flag = self._convertflag_k(self.x)
                            self.x += 1
                            self.ratio = 1.0
                
                # 处理K2面
                for i in range(-l + 1, l, 2):
                    self.n1 = i
                    for j in range(-l + 1, l, 2):
                        self.n2 = j
                        for k in range(-l + 1, l, 2):
                            self.n3 = k
                            self.r = math.sqrt(self.n1**2 + self.n2**2 + self.n3**2)
                            self.a += -1 * self.R / self.r
            else:
                # 处理K2面
                for i in range(-l, l + 1, 2):
                    self.n1 = i
                    for j in range(-l, l + 1, 2):
                        self.n2 = j
                        for k in range(-l, l + 1, 2):
                            self.n3 = k
                            self.r = math.sqrt(self.n1**2 + self.n2**2 + self.n3**2)
                            self.ratio = 1.0
                            if self.n1 == l or self.n1 == -l:
                                self.ratio /= 2
                            if self.n2 == l or self.n2 == -l:
                                self.ratio /= 2
                            if self.n3 == l or self.n3 == -l:
                                self.ratio /= 2
                            self.a += -1 * self.ratio * self.R / self.r
                            self.ratio = 1.0
                
                # 处理C60与K1面
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
                            if self.n1**2 + self.n2**2 + self.n3**2 == 0:
                                self.flag = self._convertflag_k(self.x)
                                self.x += 1
                                continue
                            self.r = math.sqrt(self.n1**2 + self.n2**2 + self.n3**2)
                            self.a += -self.flag * self.R / self.r
                            self.flag = self._convertflag_k(self.x)
                            self.x += 1
            
            self.a_array.append(self.a)
            print(f"m={self.m}\ta={self.a}")
            if self.cnt % 4 == 0:
                self.aver = sum(self.a_array) / self.cnt
                print(f"K1 result's average ={self.aver}")
            self.a = 0.0
        
        print(f"\nK1 final result's ={self.aver}\n")

    def result_c60(self, n):
        for l in range(2, n // 2 + 1):
            self.l = l
            self.m = 2 * l
            self.cnt += 1
            if l % 2 == 0:
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
                            if self.n1**2 + self.n2**2 + self.n3**2 == 0:
                                self.flag = self._convertflag_c60(self.x)
                                self.x += 1
                                continue
                            self.r = math.sqrt(self.n1**2 + self.n2**2 + self.n3**2)
                            self.ratio = 1.0
                            if self.n1 == l or self.n1 == -l:
                                self.ratio /= 2
                            if self.n2 == l or self.n2 == -l:
                                self.ratio /= 2
                            if self.n3 == l or self.n3 == -l:
                                self.ratio /= 2
                            self.a += -self.flag * self.ratio * self.R / self.r
                            self.flag = self._convertflag_c60(self.x)
                            self.x += 1
                            self.ratio = 1.0
                
                for i in range(-l + 1, l, 2):
                    self.n1 = i
                    for j in range(-l + 1, l, 2):
                        self.n2 = j
                        for k in range(-l + 1, l, 2):
                            self.n3 = k
                            self.r = math.sqrt(self.n1**2 + self.n2**2 + self.n3**2)
                            self.a += 3 * self.R / self.r
            else:
                for i in range(-l, l + 1, 2):
                    self.n1 = i
                    for j in range(-l, l + 1, 2):
                        self.n2 = j
                        for k in range(-l, l + 1, 2):
                            self.n3 = k
                            self.r = math.sqrt(self.n1**2 + self.n2**2 + self.n3**2)
                            self.ratio = 1.0
                            if self.n1 == l or self.n1 == -l:
                                self.ratio /= 2
                            if self.n2 == l or self.n2 == -l:
                                self.ratio /= 2
                            if self.n3 == l or self.n3 == -l:
                                self.ratio /= 2
                            self.a += 3 * self.ratio * self.R / self.r
                            self.ratio = 1.0
                
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
                            if self.n1**2 + self.n2**2 + self.n3**2 == 0:
                                self.flag = self._convertflag_c60(self.x)
                                self.x += 1
                                continue
                            self.r = math.sqrt(self.n1**2 + self.n2**2 + self.n3**2)
                            self.a += -self.flag * self.R / self.r
                            self.flag = self._convertflag_c60(self.x)
                            self.x += 1
            
            self.a_array.append(self.a)
            print(f"m={self.m}\ta={self.a}")
            if self.cnt % 4 == 0:
                self.aver = sum(self.a_array) / self.cnt
                print(f"C60 result's average ={self.aver}")
            self.a = 0.0
        
        print(f"\nC60 final result's ={self.aver}\n")

    def result_k2(self, n):
        for l in range(2, n // 2 + 1):
            self.l = l
            self.m = 2 * l
            self.cnt += 1
            if l % 2 == 0:
                for i in range(-l, l + 1, 2):
                    self.n1 = i
                    for j in range(-l, l + 1, 2):
                        self.n2 = j
                        for k in range(-l, l + 1, 2):
                            self.n3 = k
                            if self.n1**2 + self.n2**2 + self.n3**2 == 0:
                                continue
                            self.r = math.sqrt(self.n1**2 + self.n2**2 + self.n3**2)
                            self.ratio = 1.0
                            if self.n1 == l or self.n1 == -l:
                                self.ratio /= 2
                            if self.n2 == l or self.n2 == -l:
                                self.ratio /= 2
                            if self.n3 == l or self.n3 == -l:
                                self.ratio /= 2
                            self.a += -1 * self.ratio * self.R / self.r
                            self.ratio = 1.0
                            self.x += 1
                
                for i in range(-l + 1, l, 2):
                    self.n1 = i
                    for j in range(-l + 1, l, 2):
                        self.n2 = j
                        for k in range(-l + 1, l, 2):
                            self.n3 = k
                            self.r = math.sqrt(self.n1**2 + self.n2**2 + self.n3**2)
                            self.flag = self._convertflag_k(self.x)
                            self.x += 1
                            self.a += -self.flag * self.R / self.r
            else:
                for i in range(-l, l + 1, 2):
                    self.n1 = i
                    for j in range(-l, l + 1, 2):
                        self.n2 = j
                        for k in range(-l, l + 1, 2):
                            self.n3 = k
                            self.r = math.sqrt(self.n1**2 + self.n2**2 + self.n3**2)
                            self.ratio = 1.0
                            if self.n1 == l or self.n1 == -l:
                                self.ratio /= 2
                            if self.n2 == l or self.n2 == -l:
                                self.ratio /= 2
                            if self.n3 == l or self.n3 == -l:
                                self.ratio /= 2
                            self.flag = self._convertflag_k(self.x)
                            self.x += 1
                            self.a += -self.flag * self.ratio * self.R / self.r
                            self.ratio = 1.0
                
                for i in range(-l + 1, l, 2):
                    self.n1 = i
                    for j in range(-l + 1, l, 2):
                        self.n2 = j
                        for k in range(-l + 1, l, 2):
                            self.n3 = k
                            if self.n1**2 + self.n2**2 + self.n3**2 == 0:
                                continue
                            self.r = math.sqrt(self.n1**2 + self.n2**2 + self.n3**2)
                            self.a += -1 * self.R / self.r
            
            self.a_array.append(self.a)
            print(f"m={self.m}\ta={self.a}")
            if self.cnt % 4 == 0:
                self.aver = sum(self.a_array) / self.cnt
                print(f"K2 result's average ={self.aver}")
            self.a = 0.0
        
        print(f"\nK2 final result's ={self.aver}\n")

if __name__ == "__main__":
    k1 = K3C60()
    c60 = K3C60()
    k2 = K3C60()

    print("K1:")
    n = int(input("please enter the width of solid(n):"))
    k1.result_k1(n)

    print("C60:")
    n = int(input("please enter the width of solid(n):"))
    c60.result_c60(n)

    print("K2:")
    n = int(input("please enter the width of solid(n):"))
    k2.result_k2(n)

    input("Press Enter to exit...")