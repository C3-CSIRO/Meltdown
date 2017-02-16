from Tkinter import *
import tkFileDialog
import tkMessageBox
import tkFont
import ttk
alph="ABCDEFGHIJKLMNOP"

class PlateRunner:
	def __init__(self, master):
		self.master = master
		master.title("PlateRunner386")

		self.wells = {}
		self.wells_chem = {}
		self.scale = 40
		self.species = []

		#plate related things
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
		
		self.plate_create()		
		
		self.platerefalph.grid(row = 1, column = 0, rowspan = 8, sticky = E)
		self.platerefnum.grid(row = 0, column = 1, columnspan = 12, sticky = S)
		self.plate.grid(row = 1, column = 1, rowspan = 8, columnspan = 12, sticky = NW)

		self.empty1 = Canvas(master, width = self.scale, height = self.scale/2, bd=0, highlightthickness=0)
		self.empty1.create_rectangle(0,0,self.scale, self.scale, width=0)
		self.empty1.grid(row=9, column=1, sticky=W)

		#quadrant related things
		self.quadrant = Canvas(master, width=self.scale*3.6, height=self.scale*3.6, bd=0, highlightthickness=0)
		self.quadrant_create(384)
		self.quadrant.grid(row=10, column = 5, rowspan = 8, columnspan = 8, sticky = E)
		
		self.quadrant.xd=0
		self.quadrant.yd=0
		self.quadrant.selected=[]

		self.quadrant.bind("<Button-1>", self.quadmousedown)
		self.quadrant.bind("<Control-1>", self.quadctrl)

		#content key related things
		Label(master, text = "Enter ident:").grid(row = 10, column = 1, sticky = W)

		self.ass_input = Entry(master)
		self.ass_input.grid(row = 11, column = 1, columnspan = 4, sticky = W)

		self.assign_btn=Button(master, text="Assign selected wells", command=self.assign)
		self.assign_btn.grid(row = 12, column = 1, columnspan = 4, sticky = W)

		self.clear_btn=Button(master, text="Clear selected wells", command=self.unassign)
		self.clear_btn.grid(row = 13, column = 1, columnspan = 4, sticky = W)

		self.key_btn=Button(master, text="View/Edit Content Key", command=self.editkey)
		self.key_btn.grid(row = 10, column = 4, columnspan = 4, sticky = W)

		self.load_btn=Button(master, text="Load Content Key", command=self.load)
		self.load_btn.grid(row = 11, column = 4, columnspan = 4, sticky = W)

		self.create_btn=Button(master, text="Create Content Map", command=self.create_map)
		self.create_btn.grid(row = 15, column = 4, columnspan = 4, sticky = W)

		self.empty2=Canvas(master, width=self.scale*0.5, height=self.scale*0.5, bd=0, highlightthickness=0)
		self.empty2.grid(row=0, column=13)

		self.master.bind("<Return>", self.enter)
		self.master.bind("<Delete>", self.welldel)
		
	#create plate graphic
	def plate_create(self):
		for c in range(12):
			self.platerefnum.create_text(self.scale*(2*c+1)/2.0, self.scale/4, text=str(c*2+1)+"/"+str(c*2+2))
			for r in range(8):
				x=self.plate.create_rectangle(self.scale*c+1, self.scale*r+1, self.scale*(c+1)-1, self.scale*(r+1)-1, fill="lightblue")
				self.wells[8*c+r]=x
				for i in range(4):
					self.wells_chem[100*i+x]=""
					self.quads_chem[100*i+x]=""

		for r in range(8):
			self.platerefalph.create_text(self.scale/4, self.scale*(2*r+1)/2.0, text=alph[r*2:r*2+2])

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

	#changing well colours
	def filling(self, canvas):
		for i in range(96):
			if self.wells[i] in canvas.dragged or self.wells[i] in canvas.selected:
				canvas.itemconfig(self.wells[i], fill="red")
			else:
				lol = False
				for j in range(4):
					if self.wells_chem[100*j+i+1] != "":
						lol = True
				if lol:
					canvas.itemconfig(self.wells[i], fill="lightgreen")
				else:
					canvas.itemconfig(self.wells[i], fill="lightblue")

	def quadfilling(self, canvas, wells):
		for i in range(4):
			if self.quads[i] in canvas.selected:
				canvas.itemconfig(self.quads[i], fill="red")
			else:
				canvas.itemconfig(self.quads[i], fill="lightblue")

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
			self.quad_update(chosen)
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

	#quadrant events
	def quadmousedown(self, event):
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
		self.quadfilling(canvas, self.quads)

	def quadctrl(self, event):
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
		self.quadfilling(canvas, self.quads)
		
	def quad_update(self, chosen):
		for i in range(4):
			self.quadrant.itemconfig(2*i+2, text=self.quads_chem[i*100+chosen])

	#assign button press
	def assign(self):
		if self.plate.selected == [] or self.quadrant.selected == []:
			tkMessageBox.showwarning("Error", "No wells selected.")
		if self.plate.legend=={}:
			tkMessageBox.showwarning("Error", "No idents exist. Consider loading a content key file.")
		else:
			y=self.ass_input.get()
			if y not in self.plate.legend:
				tkMessageBox.showwarning("Error", "Ident does not exist you wet napkin!!")
			else:
				for i in self.plate.selected:
					for k in self.quadrant.selected:
						if k%2==1:
							self.quadrant.itemconfig(k+1, text=y)
							to_add = ""
							for j in self.plate.legend[y]:
								to_add+=j+"\t"
							self.wells_chem[i+50*(k-1)]=to_add[:-2]
							self.quads_chem[i+50*(k-1)]=y
		self.ass_input.delete(0,END)

	def unassign(self):
		for i in self.plate.selected:
			for k in self.quadrant.selected:
				if k%2==1:
					self.quadrant.itemconfig(k+1, text="")
					self.wells_chem[i+50*(k-1)]=""
					self.quads_chem[i+50*(k-1)]=""

	def enter(self, event):
		self.assign()

	def welldel(self, event):
		self.unassign()

	#create map
	def create_map(self):
		content_str = "Well	Condition Variable 1	Condition Variable 2	pH	d(pH)/dT	Control\n"
		for i in self.wells_chem:
			abc=alph[(i%100)%8*2-2+int(i/100)%2]
			xyz=str(int((i%100-1)/8)*2+int(i/200)+1)
			content_str+=abc
			content_str+=xyz
			content_str+="\t"
			if self.wells_chem[i]!="":
				content_str+=self.wells_chem[i]
			else:
				content_str+="EMPTY\t\t\t\t"
			content_str+="\n"
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

	#open content key window
	def load(self):
		self.plate.legend_raw=[]
		f = tkFileDialog.askopenfile(mode='r', **{'filetypes':[('all files', '.*'), ('text files', '.txt')]})
		if f != None:
			for i in f:
				self.plate.legend_raw.append(i.split("\t"))
			if self.plate.legend_raw[0][0]=="Ident":
				self.plate.legend_raw.remove(self.plate.legend_raw[0])
			for i in self.plate.legend_raw:
				self.plate.legend[i[0]]=i[1:]
			f.close()
		try:
			self.key_edit.destroy()
		except:
			None
		self.editkey()

	def editkey(self):
		try:
			self.key_edit.destroy()
		except:
			None
		self.map_header = ["Ident", "Variable 1", "Variable 2", "pH", "d(pH)/dT", "Control?"]
		self.species = []

		for i in self.plate.legend:
			self.lmao = [i]
			for j in self.plate.legend[i]:
				self.lmao.append(j)
			if self.lmao[-1][:-1]=="1":
				self.lmao[-1]="Yes"
			self.species.append(self.lmao)

		self.tree = None
		self.setup_widgets()
		self.build_tree()

		try:
			self.add_to_key.tkraise(aboveThis=None)
		except:
			None

	def setup_widgets(self):
		self.key_edit = Toplevel(self.master)
		self.key_edit.minsize(500,500)
		self.key_edit.tkraise(aboveThis=None)
		self.key_edit.title("PlateRunner Content Key")
		self.container = Frame(self.key_edit)
		self.container.grid(row=0, column=0, columnspan = 100, sticky='nesw')
		self.container.grid_columnconfigure(0, weight=5)
		self.container.grid_rowconfigure(0, weight=5)

		self.tree = ttk.Treeview(self.key_edit, columns=self.map_header, show="headings")
		self.vsb = ttk.Scrollbar(self.key_edit, orient="vertical", command=self.tree.yview)
		self.hsb = ttk.Scrollbar(self.key_edit, orient="horizontal", command=self.tree.xview)
		self.tree.configure(yscrollcommand=self.vsb.set, xscrollcommand=self.hsb.set)
		self.tree.grid(column=0, row=0, sticky='nsew', in_ = self.container)
		self.vsb.grid(column=1, row=0, sticky='ns', in_ = self.container)
		self.hsb.grid(column=0, row=1, sticky='ew', in_ = self.container)

		self.key_edit.grid_columnconfigure(0, weight=1)
		self.key_edit.grid_rowconfigure(0, weight=1)

		self.load_new_btn=Button(self.key_edit, text="Load new key", command=self.load)
		self.load_new_btn.grid(row=1,column=10,sticky=NW)

		self.container2 = Frame(self.key_edit)

		Label(self.container2, text = "Add new species to key:").grid(row = 1, column = 0, columnspan = 10, sticky = W) 

		Label(self.container2, text="Ident", anchor = W).grid(row = 2, column = 0, sticky = W)
		Label(self.container2, text="Variable 1", anchor = W).grid(row = 3, column = 0, sticky = W)
		Label(self.container2, text="Variable 2", anchor = W).grid(row = 4, column = 0, sticky = W)
		Label(self.container2, text="pH", anchor = W).grid(row = 5, column = 0, sticky = W)
		Label(self.container2, text="d(pH)/dT", anchor = W).grid(row = 6, column = 0, sticky = W)
		Label(self.container2, text="Control?",  anchor = W).grid(row = 7, column = 0, sticky = W)

		self.add_a = Entry(self.container2)
		self.add_a.grid(row=2,column=1,sticky=W)

		self.add_b = Entry(self.container2)
		self.add_b.grid(row = 3, column = 1, sticky = W)

		self.add_c = Entry(self.container2)
		self.add_c.grid(row = 4, column = 1, sticky = W)

		self.add_d = Entry(self.container2)
		self.add_d.grid(row = 5, column = 1, sticky = W)

		self.add_e = Entry(self.container2)
		self.add_e.grid(row = 6, column = 1, sticky = W)

		self.ctr = IntVar()
		self.chkbtn = Checkbutton(self.container2, variable=self.ctr)
		self.chkbtn.grid(row=7, column = 1, sticky = W)

		self.add_btn=Button(self.container2, text="Add to key", command=self.add_now_to_key)
		self.add_btn.grid(row=8,column=0,sticky=W)

		self.container2.grid(row=1, column=0, sticky = W)

	def build_tree(self):
		for col in self.map_header:
			self.tree.heading(col, text=col.title(), command=lambda c=col: self.sortby(self.tree, c, 0))
			self.tree.column(col, width=tkFont.Font().measure(col.title()))

		for item in self.species:
			self.tree.insert('', 'end', values=item)
			for ix, val in enumerate(item):
				col_w = tkFont.Font().measure(val)
				if self.tree.column(self.map_header[ix],width=None)<col_w:
					self.tree.column(self.map_header[ix], width=col_w)

	def sortby(self, tree, col, descending):
	    data = [(tree.set(child, col), child) \
	        for child in tree.get_children('')]
	    data.sort(reverse=descending)
	    for ix, item in enumerate(data):
	        tree.move(item[1], '', ix)
	    tree.heading(col, command=lambda col=col: self.sortby(tree, col, \
	        int(not descending)))

	def add_now_to_key(self):
		if self.ctr.get() == 1:
			x="1"
		else:
			x=""
		if self.add_a.get()=="":
			tkMessageBox.showinfo("Error", "Ident cannot be blank!")
			self.key_edit.lift()
		elif self.add_a.get() in self.plate.legend:
			tkMessageBox.showinfo("Overwrite?", "Ident already exists! Overwrite?")
			self.key_edit.lift()
		else:
			self.plate.legend[self.add_a.get()] = [self.add_b.get(), self.add_c.get(), self.add_d.get(), self.add_e.get(), x]
			try:
				self.key_edit.destroy()
			except:
				None
			self.editkey()

	def lol_fnc_inf(self, event):
		self.add_now_to_key()

root = Tk()

main_gui = PlateRunner(root)

root.mainloop()