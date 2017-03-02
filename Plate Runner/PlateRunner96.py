"""
PLATE RUNNER V.2


                   ***      ***                                 
                  ****      ****                                
            ------****------****------                          
           -*************************-                          
           -*************************-         *                
           -*************************-        *                 
           -*************************-        *                 
          *-*************************-***    *                  
         * -*************************-   ****                   
        *  -*************************-                          
         * -*************************-                          
         * ---------------------------                          
          *         *         **                                
           *       *            **                              
                  *              *                             
                 *              *                              
    ***         *              *                               
      ******   *              *                               
            ***              *                                
                            *                                 
                           *                                 
                          *                                  

"""

from Tkinter import *
import tkFileDialog
import tkMessageBox
import tkFont
import ttk
import re
alph="ABCDEFGHIJKLMNOP"

class PlateRunner:
	def __init__(self, master):
		self.master = master
		master.title("PlateRunner386")

		self.wells = {}
		self.wells_chem = {}
		self.scale = 50
		self.species = []

		#plate related things
		#self.mainplate = Frame(master)
		self.plate = Canvas(master, width=self.scale*12, height=self.scale*8, bd=0, highlightthickness=0)
		self.platerefalph = Canvas(master, width=self.scale/2, height=self.scale*8, bd=0, highlightthickness=0)
		self.platerefnum=Canvas(master, width=self.scale*12, height=self.scale/2, bd=0, highlightthickness=0)

		self.plate.xd=0
		self.plate.yd=0
		self.plate.selected=[]
		self.plate.dragged=[]
		self.plate.legend_raw=[]
		self.plate.legend={}

		self.quads = {}
		self.quads_chem = {}

		#event binds for plate
		self.plate.bind("<Button-1>", self.platemousedown)
		self.plate.bind("<ButtonRelease-1>", self.platemouseup)
		self.plate.bind("<B1-Motion>", self.platemousemove)
		self.plate.bind("<Control-1>", self.platectrl)
		
		#creating physical display of plate
		self.plate_create()		
		self.platerefalph.grid(row = 1, column = 0, sticky = E)
		self.platerefnum.grid(row = 0, column = 1, columnspan = 3, sticky = S)
		self.plate.grid(row = 1, column = 1, rowspan = 1, columnspan = 3, sticky = NW)

		#creating empty space
		self.empty1 = Canvas(master, width = self.scale, height = self.scale/2, bd=0, highlightthickness=0)
		self.empty1.create_rectangle(0,0,self.scale, self.scale, width=0)
		self.empty1.grid(row=3, column=1, sticky=W)

		self.empty2 = Canvas(master, width = self.scale, height = self.scale/2, bd=0, highlightthickness=0)
		self.empty2.create_rectangle(0,0,self.scale, self.scale, width=0)
		self.empty2.grid(row=6, column=1, sticky=W)


		"""
		self.quadrant = Canvas(master, width=self.scale*3.6, height=self.scale*3.6, bd=0, highlightthickness=0)
		self.quadrant_create(384)
		self.quadrant.grid(row=4, column = 3, sticky = NE)
		
		self.quadrant.xd=0
		self.quadrant.yd=0
		self.quadrant.selected=[]

		self.quadrant.bind("<Button-1>", self.quadmousedown)
		self.quadrant.bind("<Control-1>", self.quadctrl)
		"""

		#Variable 1 tree
		self.var1frame = Frame(master)

		self.var1 = ttk.Treeview(self.var1frame, columns = ["v1"], show = "headings")
		self.var1.heading("v1", text = "Variable 1")
		self.var1.column("v1", width = int(self.scale*3))

		self.var1.insert("", "end", values = ('"*No Variable 1*"'))

		self.vsb1 = ttk.Scrollbar(master, orient="vertical", command=self.var1.yview)
		self.hsb1 = ttk.Scrollbar(master, orient="horizontal", command=self.var1.xview)
		self.var1.configure(yscrollcommand=self.vsb1.set, xscrollcommand=self.hsb1.set)
		self.var1.grid(sticky=NSEW)
		self.vsb1.grid(column=1, row=0, sticky= NS, in_ = self.var1frame)
		self.hsb1.grid(column=0, row=1, sticky= EW, in_ = self.var1frame)

		self.var1_input = Entry(self.var1frame)
		self.var1_input.grid(row = 2, columnspan = 4, sticky = W)

		self.var1_btn=Button(self.var1frame, text="Add Variable 1", command=self.updatevar1)
		self.var1_btn.grid(row = 3, columnspan = 4, sticky = W)

		self.var1frame.grid(row = 4, column=1, rowspan = 2, sticky = NW)


		#Variable 2 tree
		self.var2frame = Frame(master)
		self.var2 = ttk.Treeview(self.var2frame, columns = ["v2"], show = "headings")
		self.var2.heading("v2", text = "Variable 2")
		self.var2.column("v2", width = int(self.scale*3))
		
		self.var2.insert("", "end", values = ('"*No Variable 2*"'))

		self.vsb2 = ttk.Scrollbar(master, orient="vertical", command=self.var2.yview)
		self.hsb2 = ttk.Scrollbar(master, orient="horizontal", command=self.var2.xview)
		self.var2.configure(yscrollcommand=self.vsb2.set, xscrollcommand=self.hsb2.set)
		self.var2.grid(sticky=NSEW)
		self.vsb2.grid(column=1, row=0, sticky= NS, in_ = self.var2frame)
		self.hsb2.grid(column=0, row=1, sticky= EW, in_ = self.var2frame)

		self.var2_input = Entry(self.var2frame)
		self.var2_input.grid(row = 2, columnspan = 4, sticky = W)

		self.var2_btn=Button(self.var2frame, text="Add Variable 2", command=self.updatevar2)
		self.var2_btn.grid(row = 3, columnspan = 4, sticky = W)

		self.var2frame.grid(row = 4, column = 2, rowspan = 2, sticky = N)


		#pH, dpH/dT and Control
		self.otherinputs = Frame(master)

		Label(self.otherinputs, text = "pH:").grid(row = 11, column = 9, sticky = SW)
		Label(self.otherinputs, text = "dpH/dT:  ").grid(row = 12, column = 9, sticky = W)
		Label(self.otherinputs, text = "Control?").grid(row = 13, column = 9, sticky = W)

		self.pH_input = Entry(self.otherinputs)
		self.pH_input.grid(row = 11, column = 10, columnspan = 2, sticky = SW)

		self.dpH_input = Entry(self.otherinputs)
		self.dpH_input.grid(row = 12, column = 10, columnspan = 2, sticky = W)

		self.ctr = IntVar() 
		self.chkbtn = Checkbutton(self.otherinputs, variable=self.ctr)
		self.chkbtn.grid(row=13, column = 10, sticky = W)

		self.otherinputs.grid(row = 4, column = 3, sticky = NE)


		#content key related things

		self.assign_btn=Button(master, text="Assign selected wells", command=self.assign)
		self.assign_btn.grid(row = 7, column = 1, sticky = W)

		self.clear_btn=Button(master, text="Clear selected wells", command=self.unassign)
		self.clear_btn.grid(row = 8, column = 1, sticky = W)

		self.create_btn=Button(master, text="CREATE CONTENT MAP", command=self.create_map)
		self.create_btn.grid(row = 8, column = 3, sticky = E)

		self.empty2=Canvas(master, width=self.scale*0.5, height=self.scale*0.5, bd=0, highlightthickness=0)
		self.empty2.grid(row=0, column=13)

		self.master.bind("<Return>", self.enter)
		self.master.bind("<Delete>", self.welldel)
		
	#create plate graphic

	def updatevar1(self):
		self.var1.insert("", "end", values = (re.escape(self.var1_input.get())))
		self.var1_input.delete(0,END)

	def updatevar2(self):
		self.var2.insert("", "end", values = (re.escape(self.var2_input.get())))
		self.var2_input.delete(0,END)

	def readvar(self, x):
		if x == 1:
			curItem = self.var1.focus()
			try:
				return str(self.var1.item(curItem)['values'][0])
			except:
				return ""
		elif x == 2:
			curItem = self.var2.focus()
			try:
				return str(self.var2.item(curItem)['values'][0])
			except:
				return ""

	def plate_create(self):
		for c in range(12):
			self.platerefnum.create_text(self.scale*(2*c+1)/2.0, self.scale/4, text=str(c+1))
			for r in range(8):
				x=self.plate.create_rectangle(self.scale*c+1, self.scale*r+1, self.scale*(c+1)-1, self.scale*(r+1)-1, fill="lightblue")
				z=self.plate.create_text(self.scale*(2*c+1)/2.0, self.scale*(2*r+1)/2.0, text="", fill="black")
				self.wells[8*c+r]=x
				self.wells_chem[x]=""

		for r in range(8):
			self.platerefalph.create_text(self.scale/4, self.scale*(2*r+1)/2.0, text=alph[r])

	"""
	def quadrant_create(self, wellno):
		self.quadrant.delete("all")
		if wellno==384:
			for i in [[0,0],[0,1],[1,0],[1,1]]:
				x=self.quadrant.create_oval((3*i[0]+1)/2.0*self.scale, (3*i[1]+1)/2.0*self.scale, (3*i[0]+3)/2.0*self.scale, (3*i[1]+3)/2.0*self.scale, fill="lightblue")
				z=self.quadrant.create_text((1.5*i[0]+1)*self.scale, (1.5*i[1]+1)*self.scale, text="", fill="black")
				self.quads[2*i[0]+i[1]]=x
				
		else:
			self.quadrant.create_oval(1*self.scale, 1*self.scale, 2.5*self.scale, 2.5*self.scale, fill="red")
		self.quadrant.create_rectangle(0,0,self.scale*3.5, self.scale*3.5)
	"""

	#changing well colours
	def filling(self, canvas):
		for i in range(96):
			if self.wells[i] in canvas.dragged or self.wells[i] in canvas.selected:
				canvas.itemconfig(self.wells[i], fill="red")
			else:
				canvas.itemconfig(self.wells[i], fill="lightblue")

	"""
	def quadfilling(self, canvas, wells):
		for i in range(4):
			if self.quads[i] in canvas.selected:
				canvas.itemconfig(self.quads[i], fill="red")
			else:
				canvas.itemconfig(self.quads[i], fill="lightblue")
	"""

	#mouse events
	def platemousedown(self, event):
		canvas = event.widget
		canvas.selected=[]
		canvas.xd = canvas.canvasx(event.x)
		canvas.yd = canvas.canvasy(event.y)
		raw=canvas.find_overlapping(canvas.xd-0.1, canvas.yd-0.1, canvas.xd+0.1, canvas.yd+0.1)
		if len(raw)>0:
			chosen=raw[0]
		else:
			chosen=None
		if chosen and chosen not in canvas.selected:
			canvas.selected.append(chosen)
			#self.quad_update(chosen)
		self.filling(canvas)
				
	def platemousemove(self, event):
		canvas = event.widget
		xnow = canvas.canvasx(event.x)
		ynow = canvas.canvasy(event.y)
		canvas.dragged=canvas.find_overlapping(canvas.xd, canvas.yd, xnow, ynow)
		self.filling(canvas)

	def platectrl(self, event):
		canvas = event.widget
		canvas.xd = canvas.canvasx(event.x)
		canvas.yd = canvas.canvasy(event.y)
		raw=canvas.find_overlapping(canvas.xd-0.1, canvas.yd-0.1, canvas.xd+0.1, canvas.yd+0.1)
		if len(raw)>0:
			chosen=raw[0]
		else:
			chosen=None
		if chosen and chosen not in canvas.selected:
			canvas.selected.append(chosen)
		elif chosen in canvas.selected:
			canvas.selected.remove(chosen)
		self.filling(canvas)
		
	def platemouseup(self, event):
		canvas = event.widget
		for i in canvas.dragged:
			if i not in canvas.selected:
				canvas.selected.append(i)
		self.filling(canvas)
		canvas.dragged=[]

	#assign button press
	def assign(self):
		if self.plate.selected == []:
			tkMessageBox.showwarning("Error", "No wells selected.")
		elif self.readvar(1) == "":
			tkMessageBox.showwarning("Error", "No Variable 1 selected to assign.")
		else:
			if self.ctr.get() == 1:
				x="1"
			else:
				x=""
			for k in self.plate.selected:
				if k%2==1:
					a = self.readvar(1) if self.readvar(1) != "*No Variable 1*" else ""
					aa = self.readvar(1) if self.readvar(1) != "*No Variable 1*" else "EMPTY"
					b = self.readvar(2) if self.readvar(2) != "*No Variable 2*" else ""
					self.plate.itemconfig(k+1, text = a[:6])
					self.wells_chem[k]=[aa, b, self.pH_input.get(), self.dpH_input.get(), x]

	def unassign(self):
		for k in self.plate.selected:
			if k%2==1:
				self.plate.itemconfig(k+1, text="")
				self.wells_chem[(k+1)/2]=""

	def enter(self, event):
		self.assign()

	def welldel(self, event):
		self.unassign()

	#create map
	def create_map(self):
		content_str = "Well	Condition Variable 1	Condition Variable 2	pH	d(pH)/dT	Control\n"
		for i in self.wells_chem:
			abc=alph[(i-1)%16/2]
			xyz=str(int(i/16)+1)
			content_str+=abc
			content_str+=xyz
			content_str+="\t"
			if self.wells_chem[i]!="":
				for j in self.wells_chem[i]:
					content_str+=j+"\t"
			else:
				content_str+="EMPTY\t\t\t\t\t"
			content_str = content_str[:-1]+"\n"
		gname=tkFileDialog.asksaveasfilename(title = "Save file as", filetypes = (('text files', '.txt'), ('all files', '.*')))

		if gname != "":
			if gname[-4:] != ".txt":
				gname+=".txt"
			g=open(gname, "w")
			g.write(content_str)
			g.close()
			tkMessageBox.showinfo("Success", "Content map successfully generated.")
		else:
			tkMessageBox.showwarning("Error", "Failed to save content map.")

root = Tk()

main_gui = PlateRunner(root)

root.mainloop()