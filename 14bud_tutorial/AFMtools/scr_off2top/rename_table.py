import glob
myxvgs=glob.glob('./*xvg')
mylist=['OW','MW','HW','H','H1', 'H2' ,'H3','O','O1','C','C1','C2','C3','C4','C5','N','N1'] 
for ixvg in myxvgs:
    if ixvg.startswith('./test'):
        #print ixvg
        continue
    columns=ixvg.strip('./Ala7_').strip('.xvg').split('_')
    if(mylist.index(columns[0])>mylist.index(columns[1])):
        print 'mv '+'./Ala7_'+columns[0]+'_'+columns[1]+'.xvg'+' ./Ala7_'+columns[1]+'_'+columns[0]+'.xvg'
        print columns[1]+' '+columns[0]
    else:
        print columns[0]+' '+columns[1]
