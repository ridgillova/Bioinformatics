import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import sqlite3
from PIL import Image, ImageTk
import io

# creating a cancel id variable to assign it when closing the message
cancel_id = None

# function to show info, warning and error messages (temporary)
def show_message(text, colour = 'black'):
    global cancel_id
    # if window not cleared
    if cancel_id:
        window.after_cancel(cancel_id)
    temp_message.config(text=text, fg = colour)
    # Set a time delay to clear the message
    cancel_id = window.after(10000, lambda: temp_message.config(text=''))

# warning message to confirm action
def warning_askokcancel_message(text):
    result = tk.messagebox.askokcancel(title = 'Warning', message = text)
    return result

# get all molecule names in database to add to drop-down menu
def get_molecule_names():
    conn = sqlite3.connect('molecules.db')
    cursor = conn.cursor()

    cursor.execute('SELECT Name FROM molecules;')
    molecule_names = []
    for row in cursor.fetchall():
        molecule_names.append(row[0])
    
    return molecule_names 

# finds selected molecule in database 
def find_molecule():
    selected_molecule = combobox.get()
    filters_molecule = []
    if selected_molecule == 'All molecules':
        display_data()
    elif selected_molecule != '':
        filters_molecule.append(f"Name = '{selected_molecule}'")
        display_data(filters = filters_molecule)
    else:
        show_message('No molecule selected.', 'orange')

# creates sql query and fetches data using that query
def fetch_data(filters=None):

    conn = sqlite3.connect('molecules.db')
    cursor = conn.cursor()

    sql_query = "SELECT * FROM molecules"
    if filters:
        sql_query += " WHERE " + " AND ".join(filters)

    cursor.execute(sql_query)
    data = cursor.fetchall()
    
    conn.close()
    return data

# get image data from the database to display on canvas based on selected molecule
def get_image_data_for_selected_item(name):
    conn = sqlite3.connect('molecules.db')
    c = conn.cursor()
    c.execute("SELECT Image FROM molecules WHERE Name=?", (name,))
    image_data = c.fetchone()[0]
    conn.close()
    return image_data

# display the 2D structure on canvas given the image data (binary)
def display_molecule_image(image_data):
    canvas.delete('all')
    image = Image.open(io.BytesIO(image_data))
    if (canvas.winfo_reqwidth(), canvas.winfo_reqheight()) != image.size:
        image = image.resize((canvas.winfo_reqwidth(), canvas.winfo_reqheight()), resample = Image.LANCZOS)
    photo = ImageTk.PhotoImage(image)
    canvas.create_image(0, 0, anchor=tk.NW, image=photo)
    canvas.image = photo

# show the 2D structure when a molecule is chosen from the treeview
def on_treeview_select(event):
    selected_item = tree.selection()[0]
    name = tree.item(selected_item)['values'][0]
    image_data = get_image_data_for_selected_item(name)
    display_molecule_image(image_data)
    canvas_text.delete('all')
    canvas_text.create_text(int(canvas_text.winfo_width())/2,15,text = name, font = ('Segoe UI', 11, 'bold'))
    canvas_text.create_text(int(canvas_text.winfo_width())/2,35, text = tree.item(selected_item)['values'][1], font = ('Segoe UI', 10))

# delete all filters from entry boxes
def delete_filters():
    for filter in entry_boxes:
        filter.delete(0, tk.END)

# sets everything to the 'default', clears all filters and unticks radiobuttons, as well as showing a message
def clear_all():
    delete_filters()
    checkbox_var.set(0)
    display_data()
    combobox.current(0)
    show_message('Filters cleared.', 'green')

# filter values being inserted into entry boxes for each of the radiobutton
def filter_rules():
    if checkbox_var.get() == 1:
        delete_filters()
        mw_entry_max.insert(0,'500')
        logp_entry_max.insert(0,'5')
        hba_entry_max.insert(0,'10')
        hbd_entry_max.insert(0,'5')
        display_data()
    elif checkbox_var.get() == 2:
        delete_filters()
        mw_entry_max.insert(0,'450')
        logd_entry_max.insert(0,'4')
        hba_entry_max.insert(0,'8')
        hbd_entry_max.insert(0,'4')
        rotbonds_entry_max.insert(0,'10')
        rings_entry_max.insert(0,'4')
        display_data()
    elif checkbox_var.get() == 3:
        delete_filters()
        mw_entry_max.insert(0, '500')
        logp_entry_max.insert(0,'5')
        hba_entry_max.insert(0,'10')
        hbd_entry_max.insert(0,'5')
        rotbonds_entry_max.insert(0,'10')
        far_entry_max.insert(0,'5')
        psa_entry_max.insert(0, '200')
        display_data()
    elif checkbox_var.get() == 4:
        delete_filters()
        rotbonds_entry_max.insert(0, '10')
        psa_entry_max.insert(0, '140')
        display_data()
    else:
        delete_filters()

# global variables to count all molecules being displayed and whether a molecule has been selected from the combobox
count_all = 0
molecule_selected = False
def display_data(button = 1, filters = None):
    global count_all
    global molecule_selected
    if not filters:
        filters = []
    else:
        molecule_selected = True

    try:
        if checkbox_var.get() == 3:
            filters.append(f'''(CASE WHEN MolWt <= {float(mw_entry_max.get())} THEN 1 ELSE 0 END +
                            CASE WHEN LogP <= {float(logp_entry_max.get())} THEN 1 ELSE 0 END + 
                            CASE WHEN Hbd <= {int(hbd_entry_max.get())} THEN 1 ELSE 0 END + 
                            CASE WHEN Hba <= {int(hba_entry_max.get())} THEN 1 ELSE 0 END + 
                            CASE WHEN Rotatable_bonds <= {int(rotbonds_entry_max.get())} THEN 1 ELSE 0 END + 
                            CASE WHEN PSA <= {float(psa_entry_max.get())} THEN 1 ELSE 0 END +
                            CASE WHEN FAR <= {int(far_entry_max.get())} THEN 1 ELSE 0 END)
                            >= 6''')
            for entry in list_filters:
                if entry not in [mw_entry_max, logp_entry_max, hbd_entry_max, hba_entry_max, rotbonds_entry_max, psa_entry_max, far_entry_max]:
                    if entry.get():
                        filters.append(f"{dict_filters_sql[entry]} {float(entry.get())}")
        else:
            for entry in list_filters:
                if entry.get():
                    filters.append(f"{dict_filters_sql[entry]} {float(entry.get())}")
        
        data = fetch_data(filters)
        # Clear previous data
        for row in tree.get_children():
            tree.delete(row)
        
        # Display the fetched data
            
        count = 0
        if count_all == 0:
            count_all = len(data)
        for row in data:
            count += 1
            row_new = row[1:4] + row[5:]
            tree.insert("", "end", text = row[0], values=row_new)
        if count == 0 and molecule_selected:
            show_message('The chosen molecule does not pass the selected filters.', 'orange')
        elif button == 1:
            show_message('Successfully displayed data.', 'green')

        count_message.config(text = f'{count}/{count_all} molecules displayed.') # replace 99 by variable (database count)

    except ValueError: # if cannot convert values in entry boxes to numeric variable
        show_message('Cannot apply filters. Values are not numeric.', 'red')
       
# create dictionary to save filter values with corresponding number
dict_saved_filters = {}

# insert values saves into entry boxes
def get_filters(f_number):
    filter_values = dict_saved_filters[f_number]
    delete_filters()
    for i, filter in enumerate(filter_values):
        entry_boxes[i].insert(0,filter)
    show_message(f"Filters {f_number} retrieved.", 'green')

# save values from entry boxes into dictionary with corresponding number
def save_filters(f_number):
    global dict_saved_filters
    filter_values = []
    for filter in entry_boxes:
        filter_values.append(filter.get())
    align = 'e'
    if f_number == 2:
        align = 'w'
    if f_number not in dict_saved_filters or (f_number in dict_saved_filters and warning_askokcancel_message('Do you want to overwrite previously saved filters?')):
        dict_saved_filters[f_number] = filter_values
        ttk.Button(window, text = f'Retrieve Filters {f_number}', command = lambda:get_filters(f_number)).grid(row = first_row + 6, column = first_column + f_number - 1, padx = 5, pady = 5, sticky = align)
        show_message(f"Filters saved. Retrieve using the 'Retrieve Filters {f_number}' button.", 'green')

# save figure or data (depending on extension)
def save_to_file(extension): 
    # set filetype, depending on the button pressed (extension)
    filetype = [('TSV files', '*.tsv'),('Text files', '*.txt')]
    if extension == '.png':
        filetype = ('PNG files', '*.png')
    # prompt the user to choose the file location and name
    path = filedialog.asksaveasfilename(defaultextension = extension, filetypes = [filetype, ("All files", "*.*")])
    if not path:
        show_message('No path given. Please try again.', 'orange')
        return
    try:
        if extension == '.png':
            selected_item = tree.selection()[0]
            name = tree.item(selected_item)['values'][0]
            image_data = get_image_data_for_selected_item(name)
            image = Image.open(io.BytesIO(image_data))
            image = image.resize((500,500), Image.LANCZOS)
            image.save(path, quality = 90)
            show_message(f'Figure successfuly saved in {path}.', 'green')
    # if no image selected from treeview, show message
    except IndexError:
        show_message('No image displayed. Please select a molecule to display its 2D structure.', 'red')
    # for the save data button
    if extension == '.tsv':
        with open(path, 'w') as out:
            out.write('ID\tName\tFormula\tSMILES\tMolecular Weight\tLogP\tLogD\tRotatable Bonds\tH-bond Donors\tH-bond Acceptors\tRings\tFused Aromatic Rings\tPolar Surface Area\n')
            for parent in tree.get_children():
                molecule_data = '\t'.join(map(str,tree.item(parent)['values'])) # converting values from tree into a list of strings and joining with a tab
                out.write(f'{molecule_data}\n')
        show_message(f'Data successfuly saved in {path}.', 'green')

# Global variable to track sorting order for each column
sort_order = {}
dict_column_types = {0:str, 1:str, 2:str, 3:int, 4:float, 5:float, 6:float, 7:int, 8:int, 9:int, 10:int, 11:int, 12:int, 13:float, 14:float}

# sorts molecules in tree view according to the selected column 
def sort_treeview_column(col_index):
    global sort_order
    # Get all item data as a list of tuples
    item_data = []
    for parent in tree.get_children(''):
        variable_type = dict_column_types[col_index]
        item_data.append((variable_type(tree.item(parent)['values'][col_index]), parent))
    # Ddtermine the sorting order for the current column
    if col_index not in sort_order or sort_order[col_index] == 'asc':
        item_data.sort()
        sort_order[col_index] = 'desc'
    else:
        item_data.sort(reverse=True)
        sort_order[col_index] = 'asc'
    
    # reinsert items into treeview in sorted order
    for index, (value, item) in enumerate(item_data):
        tree.move(item, '', index)


def on_configure(event):
    canvas.configure(scrollregion=canvas.bbox("all"))

#################################### WINDOW ########################################
        
window = tk.Tk()
#window.configure(bg = '#EBFCED')

window.title('Molecule Database')
window_width = 1600
window_height = 1100

screen_width = window.winfo_screenwidth()
screen_height = window.winfo_screenheight()

font_title = ('Segoe UI', 12, 'bold')
font_big = ('Segoe UI',11)
font_small = ('Segoe UI',9)

if window_height > screen_height or window_width > screen_width:
    window_width = int(screen_width * 0.9)
    window_height = int(screen_height * 0.9)
#    font_big = ('Segoe UI', int(11*0.9))
#    font_small = ('Segoe UI', int(9*0.9))
#window.geometry('1650x1100+10+20')
window.geometry(f'{window_width}x{window_height}+50+50')

first_row = 0
first_column = 0


# RULES - FILTERS

tk.Label(window, text = 'Rules', font = font_title).grid(row = first_row, column = first_column, columnspan = 2, padx = 5, pady = 5) # , sticky = 'w'
#tk.Label(window, text = '').grid(row = first_row + 1, rowspan = 6, column = first_column)
rules_values = {"Lipinski's Rules":1, "Lead-Likeness":2, "Bioavailability":3, "Veber's Rules":4} # , "GHOSE":4, "Veber's Rules":5
checkbox_var = tk.IntVar()
for key, value in rules_values.items():
    tk.Radiobutton(window, text = key, font = font_big, command = filter_rules, variable = checkbox_var, value = value).grid(row = first_row + value, column = 0, columnspan = 2, padx = 5, pady = 5, ipadx = 100) # , sticky = 'w' 

# add save filters buttons...
for i in [1,2]:
    align = 'e'
    if i == 2:
        align = 'w'
    ttk.Button(window, text = f'Save Filters {i}', command = lambda f_number = i:save_filters(f_number)).grid(row = first_row + 5, column = first_column + i - 1, padx = 5, pady = 10, ipadx = 9, sticky = align)


### FILTERS ####

tk.Label(window, text = 'Filters', font = font_title).grid(row = first_row, column = first_column + 3, columnspan = 6)

# create entry boxes (min and max) for each filter 
labels_filters = ['Molecular Weight (MW)', 'LogP', 'LogD', 'Hydrogen-bond\nDonors (HBD)', 'Hydrogen-bond\nAcceptors (HBA)', 'Number of Atoms', 'Rotatable Bonds', 'Number of Rings', ' Number of Aromatic Rings', 'Fused Aromatic Rings (FAR)', 'Polar Surface Area (PSA)', 'Quantitative Estimate\nof Drug-likeness (QED)']
entry_boxes = []
for i, label in enumerate(labels_filters):
    if i <= 5:
        row = i + first_row + 1
        column = first_column + 2
    else:
        row = i + first_row + 1 - 6
        column = first_column + 6
    ttk.Label(window, text = label, font = font_small).grid(row = row, column = column) # , padx = 5, pady = 5
    entry_min = ttk.Entry(window)
    entry_min.grid(row = row, column = column+1, padx = 5, pady = 5, sticky = 'ew')
    entry_boxes.append(entry_min)
    ttk.Label(window, text = 'to', font = font_small).grid(row = row, column = column+2) # , padx = 5, pady = 5
    entry_max = ttk.Entry(window)
    entry_max.grid(row = row, column = column+3, padx = 5, pady = 5, sticky = 'ew')
    entry_boxes.append(entry_max)

# saved Entry objects into list and in dictionary with corresponding SQL query
mw_entry_min, mw_entry_max, logp_entry_min, logp_entry_max, logd_entry_min, logd_entry_max, hbd_entry_min, hbd_entry_max, hba_entry_min, hba_entry_max, atoms_entry_min, atoms_entry_max, rotbonds_entry_min, rotbonds_entry_max, rings_entry_min, rings_entry_max, arom_rings_entry_min, arom_rings_entry_max, far_entry_min, far_entry_max, psa_entry_min, psa_entry_max, qed_entry_min, qed_entry_max = entry_boxes
list_filters = [mw_entry_min, mw_entry_max, logp_entry_min, logp_entry_max, logd_entry_min, logd_entry_max, hbd_entry_min, hbd_entry_max, hba_entry_min, hba_entry_max, atoms_entry_min, atoms_entry_max, rotbonds_entry_min, rotbonds_entry_max, rings_entry_min, rings_entry_max, arom_rings_entry_min, arom_rings_entry_max, far_entry_min, far_entry_max, psa_entry_min, psa_entry_max, qed_entry_min, qed_entry_max]
dict_filters_sql = {mw_entry_min: "MolWt >=", mw_entry_max: "MolWt <=", logp_entry_min: "LogP >=", logp_entry_max: "LogP <=", logd_entry_min: "LogD >=", logd_entry_max: "LogD <=", hbd_entry_min: "Hbd >=", hbd_entry_max: "Hbd <=", hba_entry_min: "Hba >=", hba_entry_max: "Hba <=", rotbonds_entry_min: "Rotatable_bonds >=", rotbonds_entry_max: "Rotatable_bonds <=", rings_entry_min: "Ring_count >=", rings_entry_max: "Ring_count <=", arom_rings_entry_min: "Aromatic_ring_count >=", arom_rings_entry_max: "Aromatic_ring_count <=", far_entry_min: "FAR >=", far_entry_max: "FAR <=", psa_entry_min: "PSA >=", psa_entry_max: "PSA <=", atoms_entry_min: "Atoms <=", atoms_entry_max: "Atoms <=", qed_entry_min: "QED >=", qed_entry_max: "QED <="}

# show message label
# label to show temporary messasages to inform/warn user
temp_message = tk.Label(text = '')
temp_message.grid(row = first_row + 7, column = first_column + 2, columnspan = 9, padx = (10,2), pady = 5)
# count label
count_message = tk.Label(text = '')
count_message.grid(row = first_row + 8, column = first_column + 11, padx = 5, pady = 5, sticky = 'e')

# drop down menu to select molecule from database
molecule_names = ['All molecules'] + get_molecule_names()
combobox = ttk.Combobox(window, values = molecule_names)
combobox.grid(row = first_row + 8, column = 2, columnspan = 8, sticky = 'we')
ttk.Button(window, text = 'Find Molecule', command = find_molecule).grid(row = first_row + 8, column = first_column + 10, padx = 5, pady = 5)


#### TREEVIEW ####

tree = ttk.Treeview(window, columns = ("Name", "Formula", "SMILES", "# Atoms","MW", "LogP", "LogD", "# RotBonds", "# HBD", "# HBA", "# Rings", "# AromRings","FAR", "PSA", "QED"))
tree.heading("#0", text = "ID")
tree.column("#0", width = 10)
columns_widths = [50, 50, 200, 30, 30, 20, 20, 50, 20, 20, 20, 50, 20, 20, 20]
for index, col_name in enumerate(tree['columns']):
    tree.heading(index, text=col_name, command = lambda col_index = index:sort_treeview_column(col_index))
    tree.column(index, width = columns_widths[index])
# adding the command when treeview item selected
tree.bind("<<TreeviewSelect>>", on_treeview_select)
tree.grid(row = first_row + 9, rowspan = 10, column = first_column + 2, columnspan = 10, pady = 5, sticky = 'nsew')
# adding a scrollbar to the treeview
scrollbar = ttk.Scrollbar(window, orient=tk.VERTICAL, command=tree.yview)
tree.configure(yscroll=scrollbar.set)
scrollbar.grid(row = first_row + 9, rowspan = 10, column = first_column + 12, sticky = 'nsw')

ttk.Button(window, text = "Apply Filters", command = lambda: display_data(1)).grid(row = first_row + 3, column = first_column + 10, columnspan = 2, padx = 5, ipadx = 10, ipady = 5, pady = 5) # , sticky = 'nsew'
ttk.Button(window, text = "Clear Filters", command = clear_all).grid(row = first_row + 4, column = first_column + 10, columnspan = 2, padx = 5, pady = 5, sticky = 's')

# create canvas for name, formula and for 2D structure
ttk.Label(window, text = 'Click on the molecule to display 2D structure:').grid(row = first_row + 8, column = 0, columnspan = 2,padx = 5, pady = 5, sticky = 'n')
canvas_text = tk.Canvas(window, width = 300, height = 50)
canvas_text.grid(row = first_row + 9, column = 0, columnspan = 2, sticky = 'nsew')
canvas = tk.Canvas(window, width = 400, height = 400)
canvas.grid(row = first_row + 10, rowspan = 4, column = 0, columnspan = 2, padx = 10, pady = 5, sticky = 'nsew')

# display cartoon molecule when opening the GUI (gives user an idea that's where the 2D structure should appear)
molecule_image = Image.open('molecule.png')
if (canvas.winfo_reqwidth(), canvas.winfo_reqheight()) != molecule_image.size:
    molecule_image = molecule_image.resize((canvas.winfo_reqwidth(), canvas.winfo_reqheight()), resample = Image.LANCZOS)
molecule_photo = ImageTk.PhotoImage(molecule_image)
molecule_image_item = canvas.create_image(0,0, anchor = 'nw', image = molecule_photo)


# buttons for saving images and data
ttk.Button(window, text = 'Save figure', command = lambda:save_to_file('.png')).grid(row = first_row + 16, column = first_column, columnspan = 2, padx = 5, pady = 20, ipadx = 5, ipady = 2) # , sticky = 's'
ttk.Button(window, text = 'Save data', command = lambda:save_to_file('.tsv')). grid(row = first_row + 17, column = first_column, columnspan = 2 , padx = 5, pady = 20, ipadx = 5, ipady = 2)

window.rowconfigure(first_row + 18, weight = 1) # makes the treeview fill up to the end of the window

# for resizing the widgets (not fully functional, does not make the text smaller and makes the rows and columns the same ratio and not according to the label length etc.)
'''
for i in range(0,19):
#    window.grid_rowconfigure(i, weight=1)
    window.rowconfigure(i, weight = 1, minsize = 35)
for i in range(0,12):
#    window.grid_columnconfigure(i, weight=1)
    window.columnconfigure(i, weight = 1, minsize = 100)
'''

# displays all molecules from the database when the GUI is opened (would be deleted with more molecules in the database - could fill up memory space)
display_data(0)

# count all molecules displayed (in database)
count_all = len(tree.get_children(''))

window.mainloop()
