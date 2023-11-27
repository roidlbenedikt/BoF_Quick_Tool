# -*- coding: utf-8 -*-
"""
Created on Sep 17  2023
Verion: 0.01
@author: 212707803
"""


import cv2
import numpy as np
import os
import csv
from matplotlib import pyplot as plt

import tkinter as tk
from tkinter import filedialog
import os
import subprocess

# Function to open a folder dialog and store the selected folder location
def open_folder():
    folder_path = filedialog.askdirectory()
    folder_path_label.config(text="Selected Folder: " + folder_path)
    run_button.config(state=tk.NORMAL)
    # Store the folder path for later use
    global selected_folder
    selected_folder = folder_path

# Function to run a program using the selected folder location and show status updates
def run_program():
    if selected_folder:
        status_text.config(state=tk.NORMAL)
        status_text.delete(1.0, tk.END)  # Clear previous status updates
        status_text.insert(tk.END, "Running program...\n")
        status_text.update()  # Update the text widget to display the message

        # Replace 'my_function' with the actual function you want to run
        program_output = check_misalignmnets(selected_folder)

        status_text.insert(tk.END, program_output)

        status_text.config(state=tk.DISABLED)


def remove_duplicates(array_tobe_cleaned):
    for keys in array_tobe_cleaned:
        relevant_codes=list(set(array_tobe_cleaned[keys]['code']))
        for codes_relevant in relevant_codes:

            index_codes = [ii for ii, x in enumerate(array_tobe_cleaned[keys]['code']) if codes_relevant == x]
            x_array= np.array(array_tobe_cleaned[keys]['x'])
            array_tobe_cleaned[keys]['x'][index_codes[0]]= np.mean(x_array[index_codes])
            y_array = np.array(array_tobe_cleaned[keys]['y'])
            array_tobe_cleaned[keys]['y'][index_codes[0]] = np.mean(y_array[index_codes])
            radial_array = np.array(array_tobe_cleaned[keys]['radial'])
            array_tobe_cleaned[keys]['radial'][index_codes[0]] = np.mean(radial_array[index_codes])

            for keys2 in array_tobe_cleaned[keys]:
                try:
                    for elem in sorted(index_codes[1:], reverse=True):
                        del(array_tobe_cleaned[keys][keys2][elem])
                except:
                    pass

    return array_tobe_cleaned


def  get_misalignment(rel_misalignment,all_radius,allcx,allcy,ref_radii,rel_radius_index, centersx, rel_fields,circle_info,conversion,ii,rad_index1,rad_index2,ref_radindex,rel_os):

    if all_radius > ref_radii[ref_radindex] * 0.9 and all_radius < ref_radii[ref_radindex] * 1.1:
             
        try:
            deviationx = round(allcx[rel_radius_index[rad_index1][rad_index2]] - centersx,4) #np.mean(allcx)
            deviationy = round(allcy[rel_radius_index[rad_index1][rad_index2]] - allcy[ii],4) #np.mean(allcy)
            radial_deviation = round(np.sqrt((deviationx * conversion) ** 2 + (deviationy * conversion) ** 2),4)
            rel_misalignment[rel_os][rel_fields[4]].append(radial_deviation)
            rel_misalignment[rel_os][rel_fields[0]].append(-deviationx*conversion)
            rel_misalignment[rel_os][rel_fields[1]].append(deviationy * conversion)
                        
            rel_misalignment[rel_os][rel_fields[5]].append(circle_info[rel_fields[5]])
            rel_misalignment[rel_os][rel_fields[6]].append(circle_info[rel_fields[6]])
            rel_misalignment[rel_os][rel_fields[7]].append(circle_info[rel_fields[7]])
            rel_misalignment[rel_os][rel_fields[8]].append(rel_os)
        except:
            rel_misalignment = rel_misalignment

        return rel_misalignment

def getmaxmin_value(mismatch_values):
    max_values=0
    max_index=0
    for i in range(len(mismatch_values)):
        if np.absolute(mismatch_values[i]) > max_values:
            max_values = np.absolute(mismatch_values[i])
            max_index = i
    return np.sign(mismatch_values[max_index])*max_values

def rename_files(path_to_folder,files_to_rename):
    files_renamed=[]
    for files in files_to_rename:
        files_replaced=files.replace(u'\u200b','')
        #files_renamed.append(files_replaced)
        os.rename(os.path.join(path_to_folder, files), os.path.join(path_to_folder, files_replaced))
    return

def get_global_positions(size_img,max_radius):

    global_positions = {}

    global_center_x = max_radius['x']#size_img[1] / 2
    global_center_y = max_radius['y']#size_img[0] / 2
    factor_high = 2
    factor_low = 0.9

    global_positions['global_center_x'] = [global_center_x * 1.3, global_center_x * 0.7]
    global_positions['global_center_y'] = [global_center_y * 1.3, global_center_y * 0.7]

    global_positions['global_upper_right_OS2_x'] = [global_center_x+max_radius['radius']*factor_low, global_center_x +max_radius['radius']*factor_high]#[size_img[1] * 0.65, size_img[1]]
    global_positions['global_upper_right_OS2_y'] = [global_center_y -max_radius['radius']*factor_high, global_center_y-max_radius['radius']*factor_low]# [0,size_img[0] * 0.25]
    global_positions['global_upper_left_OS1_x'] = [global_center_x -max_radius['radius']*factor_high, global_center_x-max_radius['radius']*factor_low]#[0,size_img[1] * 0.25]
    global_positions['global_upper_left_OS1_y'] = [global_center_y -max_radius['radius']*factor_high, global_center_y-max_radius['radius']*factor_low]#[0,size_img[0] * 0.25]

    global_positions['global_lower_right_OS4_x'] = [global_center_x+max_radius['radius']*factor_low, global_center_x +max_radius['radius']*factor_high]#[size_img[1] * 0.65, size_img[1]]
    global_positions['global_lower_right_OS4_y'] = [global_center_y+max_radius['radius']*factor_low, global_center_y +max_radius['radius']*factor_high]#[size_img[0] * 0.65, size_img[0]]
    global_positions['global_lower_left_OS3_x'] =  [global_center_x - max_radius['radius']*factor_high, global_center_x - max_radius['radius']*factor_low]#[0, size_img[1] * 0.25]
    global_positions['global_lower_left_OS3_y'] =  [global_center_y + max_radius['radius']*factor_low, global_center_y + max_radius['radius']*factor_high]#[size_img[0] * 0.65, size_img[0]]

    return global_positions

def check_misalignmnets(path_to_folder):

    #input
    #path_to_folder = r"C://Users//212707803//Box//05 Mline//LPS14//2023_10_04//Bilder//"
    #r"C://Users//212707803//Desktop//relevant_GE//Technology//Robust Stitching//Stitching_General//SubProjects//WP1_Misalignment_2020//Data//fulldata//Process_raw_data//Burn-on-foil//Bilder//"#
    #r"C://Users//212707803//Box//MLINE - Aviation Additive Collaboration//08_Stitching//02_Machine_capability//02_Documentation//Post-processing routines//For Burn on Foil//GEADDL-MLINELPS-900017_BoF after power outage of 30-09-2023//"#C:/Users/212707803/Downloads/BoF/"
    #r"C://Users//212707803//Downloads//BoF//BoF_smartphone//test//"
    #r"C:/Users/212707803/Downloads/BoF/BoF/"#"./"#r"C:/Users/212707803/Downloads/BoF/"
    #r"C://Users//212707803//Desktop//relevant_GE//Technology//Robust Stitching//Stitching_General//SubProjects//WP1_Misalignment_2020//Data//fulldata//Process_raw_data//Burn-on-foil//Bilder//"#

    output_path=path_to_folder+"//output"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    ref_radii=[3,2.5,2,1.5]#[3+0*0.025,2.5+0*0.025,2+0*0.025,1.5+0*0.025]# # outer to inner, 3mm most outer cirlce diameter, OS4, OS3, OS2, OS1
    pic_count=1
    blackwhite_thresholds=[160, 255]
    blur_parameters=[11, 11] # in Pixel
    misalignment_limit = 0.5 # in mm
    factor_to_look=10 #which fraction of the max. radius is interesing -> corner findings
    rel_combinations=['OS1-OS3','OS2-OS4','OS1-OS2','OS2-OS3','OS3-OS4','OS1-OS4']
    output_csv="Output.csv"



    #write input on screen
    print("----------------------------------")
    print("relevant radii in mm: ",*ref_radii)
    print("Black&White threshold: ",*blackwhite_thresholds)
    print("Blur parameters: ",*blur_parameters)
    print("Max. allowed mismatch in mm "+str(misalignment_limit))
    print("relevant OS-to-OS combinations: ",*rel_combinations)
    print("----------------------------------")

    rel_misalignment={}
    rel_fields=['x','y','xmax','ymax','radial','posx','posy','code','OS']
    failed_pictures=[]

    #initializing rel_misalignment
    for keys in rel_combinations:
        rel_misalignment[keys]={}
        rel_misalignment[keys][rel_fields[0]] = []
        rel_misalignment[keys][rel_fields[1]] = []
        rel_misalignment[keys][rel_fields[2]] = []
        rel_misalignment[keys][rel_fields[3]] = []
        rel_misalignment[keys][rel_fields[4]] = []
        rel_misalignment[keys][rel_fields[5]] = []
        rel_misalignment[keys][rel_fields[6]] = []
        rel_misalignment[keys][rel_fields[7]] = []
        rel_misalignment[keys][rel_fields[8]] = []

    onlyfiles = [f for f in os.listdir(path_to_folder) if os.path.isfile(os.path.join(path_to_folder, f))]
    #output_csv.append([row], [column], [x], [y], )

#   Dynolite might throw fancy \u200b spaces
    rename_files(path_to_folder,onlyfiles)

    for filename in onlyfiles:#os.listdir(path_to_folder):
        if filename.endswith((".jpg",".png",".tiff",".tif",".JPG",".PNG",".TIFF",".TIF")):
            img = cv2.imread(os.path.join(path_to_folder, filename)) #
            print("processing picture: " + os.path.join(path_to_folder, filename))

            # get code, position, etc. from the file
            # OS_involved=filename.split("_")
            circle_info = {}


            circle_info['all'] = filename.replace('.', '_').split('_')
            if filename[0:2] == "GE":
                circle_info[rel_fields[7]] = circle_info['all'][1] + '_' + circle_info['all'][2]
            else:
                circle_info[rel_fields[7]] = circle_info['all'][0]+ '_' + circle_info['all'][1]#circle_info['all'][0] + '_' + circle_info['all'][1]
            try:
                circle_info[rel_fields[5]] = circle_info['all'][2].split('x')[1]
                circle_info[rel_fields[6]] = circle_info['all'][3].split('y')[1]
            except:
                circle_info[rel_fields[5]] = 0
                circle_info[rel_fields[6]] = 0

            # convert the image to grayscale
            gray = cv2.cvtColor(np.array(img), cv2.COLOR_RGB2GRAY)
            size_img = np.shape(np.array(img))
        # get center region and corners:

            OS_found={}
            OS_found['OS1']=0
            OS_found['OS2'] = 0
            OS_found['OS3'] = 0
            OS_found['OS4'] = 0
            OS_found['Total']=0


        #Assign radiii
            rel_radii=[]
           # for OS in OS_involved:
           #     if OS == "OS4":
                    #rel_radii.append(ref_radii[0])
           #     elif OS == "OS3":
                    #rel_radii.append(ref_radii[1])
           #     elif OS == "OS2":
                    #rel_radii.append(ref_radii[2])
           #     elif OS == "OS1":
                    #rel_radii.append(ref_radii[3])

            blur = cv2.GaussianBlur(gray, (blur_parameters[0], blur_parameters[1]), 0)

            thresh = np.uint8(cv2.threshold(blur, blackwhite_thresholds[0], blackwhite_thresholds[1] ,cv2.THRESH_BINARY+cv2.THRESH_OTSU)[1])
            canny =cv2.Canny(thresh, 0, 1)

            contours, hierarchy = cv2.findContours(canny,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)

            closed_contours = []
            open_contours = []
            for i in contours:
                if cv2.contourArea(i) > cv2.arcLength(i, True):
                    closed_contours.append(i)

            #go forward with closed contours only
            contours = closed_contours
            #cv2.drawContours(img, closed_contours, -1, (0,255,0), 1)
            #cv2.imshow("Image", img)
            #cv2.waitKey()
            #plt.imshow(gray_image)

            allareaC = []
            allareaV = []
            allcx=[]
            allcy=[]
            all_cradius=[]

            for ii,c in enumerate(contours):
                # calculate moments for each contour
                (x, y), radius = cv2.minEnclosingCircle(c)
                center = (int(x), int(y))
                radius = int(radius)
                all_cradius.append(radius)

            max_radius={}
            max_radius['radius']=np.max(np.array(all_cradius))
            max_radius_index=np.argmax(np.array(all_cradius))

            M = cv2.moments(contours[max_radius_index])
            # calculate x,y coordinate of center
            cX = int(M["m10"] / M["m00"])
            cY = int(M["m01"] / M["m00"])
            max_radius['x'] = cX
            max_radius['y'] = cY

            global_positions = get_global_positions(size_img, max_radius)

            all_radius = []
            OS_found['relevantOS']=[]
            for c_count,c in enumerate(contours):
                # calculate moments for each contour
                #(x, y), radius = cv2.minEnclosingCircle(c)
                #center = (int(x), int(y))
                radius = int(all_cradius[c_count])
                if radius > max_radius['radius']/factor_to_look: # make that smaller to find the corners
                    allareaC.append(float(radius) ** 2 * np.pi)
                    allareaV.append(cv2.contourArea(c))

                    M = cv2.moments(c)
                    # calculate x,y coordinate of center
                    try:
                        cX = int(M["m10"] / M["m00"])
                        cY = int(M["m01"] / M["m00"])
                    except:
                        cX=0
                        cY=0

                    #cv2.circle(img, (cX, cY), 5, (255, 255, 255), -1)
                    #cv2.drawContours(img, c, -1, (0,255,0), 3)
                    #cv2.imshow("Image", img)
                    #cv2.waitKey()

                    # get L s in the corner
                    #OS1
                    if cX < global_positions['global_upper_left_OS1_x'][1] and cX > global_positions['global_upper_left_OS1_x'][0] and \
                            cY < global_positions['global_upper_left_OS1_y'][1] and cY > global_positions['global_upper_left_OS1_y'][0] and OS_found['OS1']==0:
                        rel_radii.append(ref_radii[3])
                        OS_found['OS1']=1
                        OS_found['Total'] += 1
                        OS_found['relevantOS'].append('OS1')
                    elif cX < global_positions['global_upper_right_OS2_x'][1] and cX > global_positions['global_upper_right_OS2_x'][0] and \
                            cY < global_positions['global_upper_right_OS2_y'][1] and cY > global_positions['global_upper_right_OS2_y'][0] and OS_found['OS2']==0:
                        rel_radii.append(ref_radii[2])
                        OS_found['OS2'] = 1
                        OS_found['Total'] += 1
                        OS_found['relevantOS'].append('OS2')
                    elif cX < global_positions['global_lower_left_OS3_x'][1] and cX > global_positions['global_lower_left_OS3_x'][0] and \
                            cY < global_positions['global_lower_left_OS3_y'][1] and cY > global_positions['global_lower_left_OS3_y'][0] and OS_found['OS3']==0:
                        rel_radii.append(ref_radii[1])
                        OS_found['OS3'] = 1
                        OS_found['Total'] += 1
                        OS_found['relevantOS'].append('OS3')
                    elif cX < global_positions['global_lower_right_OS4_x'][1] and cX > global_positions['global_lower_right_OS4_x'][0] and \
                            cY < global_positions['global_lower_right_OS4_y'][1] and cY > global_positions['global_lower_right_OS4_y'][0] and OS_found['OS4']==0:
                        rel_radii.append(ref_radii[0])
                        OS_found['OS4'] = 1
                        OS_found['Total'] += 1
                        OS_found['relevantOS'].append('OS4')


                    # circle centers are clustered in the center of the picture
                    if cX < global_positions['global_center_x'][0] and cX > global_positions['global_center_x'][1] and \
                            cY < global_positions['global_center_y'][0] and cY > global_positions['global_center_y'][1]:

                        allcx.append(cX)
                        allcy.append(cY)
                        all_radius.append(radius)

            if all_radius and rel_radii:

                conversion = max(rel_radii) / max(all_radius)  # assuming that the most outer circle is compared to the right diameter according to the first OS in the name
                # get an idea about the reference center
                ref_centerx = np.median(allcx)
                ref_centery = np.median(allcy)

                rel_radius_index=[]
                for ii in range(0,4):
                    rel_radius_index.append([])

                for ii,centersx in enumerate(allcx):
                    # look for ALL circles
                    if all_radius[ii] * conversion > ref_radii[0] * 0.9 and all_radius[ii] * conversion < ref_radii[0] * 1.1:
                        rel_radius_index[0].append(ii)
                    elif all_radius[ii] * conversion > ref_radii[1] * 0.9 and all_radius[ii] * conversion < ref_radii[1] * 1.1:
                        rel_radius_index[1].append(ii)
                    elif all_radius[ii] * conversion > ref_radii[2] * 0.9 and all_radius[ii] * conversion < ref_radii[2] * 1.1:
                        rel_radius_index[2].append(ii)
                    elif all_radius[ii] * conversion > ref_radii[3] * 0.9 and all_radius[ii] * conversion < ref_radii[3] * 1.1:
                        rel_radius_index[3].append(ii)



                for ii,centersx in enumerate(allcx):

                    #Only sensical misalignmnents
                    if np.sqrt(np.abs(centersx - np.mean(allcx))**2+np.abs(allcy[ii] - np.mean(allcy))**2)*conversion<misalignment_limit and OS_found['Total'] >1:
                        rel_os = ""
                        get_misalignment(rel_misalignment,all_radius[ii]*conversion,allcx,allcy,ref_radii,rel_radius_index, centersx,rel_fields, circle_info,conversion,ii,0,0,1,'OS3-OS4')
                        get_misalignment(rel_misalignment,all_radius[ii]*conversion,allcx,allcy,ref_radii,rel_radius_index, centersx,rel_fields, circle_info,conversion,ii,1,0,2,'OS2-OS3')
                        get_misalignment(rel_misalignment,all_radius[ii]*conversion,allcx,allcy,ref_radii,rel_radius_index, centersx,rel_fields, circle_info,conversion,ii,1,0,3,'OS1-OS3')
                        get_misalignment(rel_misalignment,all_radius[ii]*conversion,allcx,allcy,ref_radii,rel_radius_index, centersx,rel_fields, circle_info,conversion,ii,0,0,2,'OS2-OS4')
                        get_misalignment(rel_misalignment,all_radius[ii]*conversion,allcx,allcy,ref_radii,rel_radius_index, centersx,rel_fields, circle_info,conversion,ii,2,0,3,'OS1-OS2')
                        get_misalignment(rel_misalignment,all_radius[ii]*conversion,allcx,allcy,ref_radii,rel_radius_index, centersx,rel_fields, circle_info,conversion,ii,0,0,3,'OS1-OS4')

            else:
                failed_pictures.append(filename)


            pic_count=pic_count+1


            print(OS_found['Total'])

    rel_misalignment = remove_duplicates(rel_misalignment)

    # display the image
    fig, axs = plt.subplots(3,2)
    axs = axs.flatten()
    fig_count=0
    for OS_combinations in rel_misalignment:
        axs[fig_count].hist(rel_misalignment[OS_combinations]['radial'], 3)
        axs[fig_count].set_title(OS_combinations)
        axs[fig_count].set_xlim([-0.2, 0.2])
        fig_count+=1
    for ax in axs.flat:
        ax.set(xlabel='Radial mismatch (mm)', ylabel='Occurences')
    fig.tight_layout()
    fig.savefig(output_path+'//Radialmismatch.png', bbox_inches="tight")

    fig, axs = plt.subplots(3,2)
    axs = axs.flatten()
    fig_count=0
    for OS_combinations in rel_misalignment:
        axs[fig_count].hist(rel_misalignment[OS_combinations]['x'], 3)
        axs[fig_count].set_title(OS_combinations)
        axs[fig_count].set_xlim([-0.2, 0.2])
        fig_count+=1
    for ax in axs.flat:
        ax.set(xlabel='X mismatch (mm)', ylabel='Occurences')
    fig.tight_layout()
    fig.savefig(output_path+'//Xmismatch.png', bbox_inches="tight")

    fig, axs = plt.subplots(3,2)
    axs = axs.flatten()
    fig_count=0
    for OS_combinations in rel_misalignment:
        axs[fig_count].hist(rel_misalignment[OS_combinations]['y'], 3)
        axs[fig_count].set_title(OS_combinations)
        axs[fig_count].set_xlim([-0.2, 0.2])
        fig_count+=1
    for ax in axs.flat:
        ax.set(xlabel='Y mismatch (mm)', ylabel='Occurences')
    fig.tight_layout()
    fig.savefig(output_path+'//Ymismatch.png', bbox_inches="tight")

    output_string=["Output written to: "+ output_path +"\n\n Maximum mismatch in (mm):\n\n"]

    # write output for csv
    for OS_combinations in rel_misalignment:
        try:
            output_string.append(str(OS_combinations+": \t x: " + "{:10.4f}".format(getmaxmin_value(rel_misalignment[OS_combinations]['x'])) + "\t y: " + "{:10.4f}".format(getmaxmin_value(rel_misalignment[OS_combinations]['y'])) + "\t radial: " + "{:10.4f}".format(max(rel_misalignment[OS_combinations]['radial'])) + "\n"))
        except:
            output_string.append(OS_combinations+": not existing \n")

    for failed in failed_pictures:
        output_string.append("Failed to analyze: "+failed+"\n")


    # text output
    print("------------------------------------")
    #print("Maximum misalignment in ")
    [print(a) for a in output_string]

    with open(output_path+"//Output.txt", "w") as text_file:
        print("------------------------------------", file=text_file)
       # print("Maximum misalignment in ", file=text_file)
        [print(a, file=text_file) for a in output_string]

    #csv-output
    fields=['x','y','radial','posx','posy','code','OS']
    with open(output_path+"//"+output_csv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        #csvwriter.writerows(rows)
        for key1 in rel_misalignment:

            for ii in range(0,len(rel_misalignment[key1]['x'])):
                output_string_csv = []

                for key2 in fields:
                    #print(rel_misalignment[key1][key2][ii],key2)
                    output_string_csv.append(rel_misalignment[key1][key2][ii])
                csvwriter.writerow(output_string_csv)


    print("CSV-output written")
    plt.show()
    return output_string


# Create the main window
root = tk.Tk()
root.title("Burn-on-foil-picture processor, V1.0")

# Create a label to display the selected folder
folder_path_label = tk.Label(root, text="Pictures location: ")
folder_path_label.pack()

# Create a button to open a folder dialog
open_button = tk.Button(root, text="Open Folder", command=open_folder)
open_button.pack()

# Create a button to run the program using the selected folder
run_button = tk.Button(root, text="Run Analysis", command=run_program, state=tk.DISABLED)
run_button.pack()

# Create a text widget to display status updates
status_text = tk.Text(root, wrap=tk.WORD, height=10, width=120)
status_text.pack()
status_text.config(state=tk.DISABLED)  # Make the text widget read-only

# Initialize the selected_folder variable
selected_folder = None

# Start the Tkinter main loop
root.mainloop()


