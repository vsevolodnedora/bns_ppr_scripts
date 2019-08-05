# from tkinter import *
from _tkinter import *
from Tkinter import Label, Button, Entry, Tk, BooleanVar, Radiobutton, IntVar, INSERT, END
import ScrolledText as tkst
import tkMessageBox
from ttk import Combobox, Checkbutton

window = Tk()

''' titile '''
window.title("Welcome to LikeGeeks app")

''' label '''
lbl = Label(window, text="Hello", font=("Arial Bold", 50))
lbl.grid(column=0, row=0) # necessary to set the location of the

''' window size '''
window.geometry('450x500')

''' button '''
def clicked():
    ''' execute this fucntion when button is clicked '''
    lbl2 = Label(window, text="Clicked!", font=("Arial Bold", 12)) # wll show the new text
    print("button 1")
    lbl2.grid(column=0, row=1)  # necessary to set the location of the
    # lbl2.configure(text="Button was clicked !!")

btn = Button(window, text="Click Me", bg="orange", fg="red", command=clicked) # with colros annd fucntion
btn.grid(column=1, row=0) # set the location


''' imput text '''
txt = Entry(window, width=10)
txt.grid(column=0, row=2)
input_text = txt.get()
txt.focus() # will focus the user on this line
# txt = Entry(window,width=10, state='disabled') # does not allow user to enter the text
# use the new button to show the text
def show_text():
    ''' show this text when the button is clicked'''
    lbl3 = Label(window, text="input text", font=("Arial", 12))
    lbl3.grid(column=2, row=2)
btn2 = Button(window, text = "show text", bg="blue", fg="red", command=show_text)
btn2.grid(column=1, row=2)

''' combox '''

combo = Combobox(window)
combo['values']= (1, 2, 3, 4, 5, "Text")
combo.current(1) #set the selected item
combo.grid(column=0, row=3)
val = combo.get()
# print(val)
chk_state = BooleanVar()
chk_state.set(False) #set check state
chk = Checkbutton(window, text='Choose', var=chk_state)
chk.grid(column=1, row=3)

''' radiobutton ( [1] [2] [3] button )'''

selected = IntVar()
rad1 = Radiobutton(window,text='First', value=1, variable=selected)
rad2 = Radiobutton(window,text='Second', value=2, variable=selected)
rad3 = Radiobutton(window,text='Third', value=3, variable=selected)
def clicked():
   print(selected.get())
btn = Button(window, text="Click Me", command=clicked)
rad1.grid(column=0, row=4)
rad2.grid(column=1, row=4)
rad3.grid(column=2, row=4)
btn.grid(column=3, row=4)

''' textarea '''

txt = tkst.ScrolledText(window, width=40, height=10)
txt.grid(column=0, row=5)
txt.insert(INSERT,'You text goes here')
def delete_text():
    ''' show this text when the button is clicked'''
    txt.delete(1.0,END) # to clear the text
btn_to_delete = Button(window, text="clear", font=("Arial", 12), command=delete_text)
btn_to_delete.grid(column=1, row=5)

''' messegebox '''
def show_messege_conent():
    tkMessageBox.showinfo('Message title', 'Message content')
btn_for_messege = Button(window, text="show window", font=("Arial", 12), command=show_messege_conent)
btn_for_messege.grid(column=0, row=6)


''' mail loop'''
window.mainloop()