import tkinter as tk

class Window:
    def __init__(self,master):
        self.master = root
        
        self.frame = tk.Frame(self.master, width = 300, height = 200)
        self.frame.pack()
        
        self.scale = tk.Scale(self.frame, from=0, to=10, orient=tk.HORIZONTAL, command=self.ret,
                              resolution=1, digit=2, tickinterval=1, length=200, sliderlength=20,
                              label="My scale", showvalue=0, troughcolor=s"blue")
        
        self.scale.place(x=30,y=30)
        
    def ret(self,value):
        print(value)
        
        
root = tk.Tk()
window = Window(root)
root.mainloop()