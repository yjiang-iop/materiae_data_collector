#!/usr/bin/env python
import os 
import numpy as np
from functools import reduce
import re

Z4gid = [52,56,58,60,61,62,70,88,126,130,133,135,136, \
       137,138,141,142,163,165,167,202,203,205,222, \
       223,227,228,230 ]

def read_indicator_table(gid, ind):
    #input: 
    #space group id: gid, indicator list: ind, in list form, e.g ind = [1,1,1,2]
    #return:
    #ind_table: lines of the indicator table which correspond to the input indicator
    #table_head: head of the table, explain the content of each column

    ind_str = str(reduce(lambda x,y: str(x)+str(y), ind)) #add str(), since if len(ind)=1, reduce will return a number, not a str
    dir = os.path.split(os.path.realpath(__file__))[0] 
    with open(dir + '/tabind.tex','r') as file:
        filelines = file.readlines()
        for cnt, line in enumerate(filelines):
            if 'Space group' in line:
                current_gid = int(line.split('#')[1].split()[0])
                if current_gid == gid:
                    count = cnt + 3
                    head = filelines[count].strip('\n\\').split('&')
                    head = [i.strip() for i in head[1:]] # remve extra space

                    count += 3
                    body = []
                    while '&' in filelines[count]: 
                        tmp = filelines[count].strip('\n\\').split('&')
                        tmp = [i.strip() for i in tmp] # remove extra space
                        
                        # for sg175,191, whose ind class are Z6*Z12, need to remove comma in the indicator read form tex file
                        if ',' in tmp[0]:
                            tmp_ind = [i.strip() for i in tmp[0].split(',')]
                            tmp[0] = reduce(lambda x,y: str(x)+str(y), tmp_ind)

                        tr = []
                        for item in tmp[1:]:
                            observable = False
                            p = r'\\textcolor\{blue\}\{\$\\mathds\{(.*)\}\$\}'
                            if re.match(p, item):
                                observable = True
                                item = re.match(p, item).group(1)
                            p = r'\$(.*)\$'
                            if re.match(p, item):
                                item = re.match(p, item).group(1) 
                            tr.append({'observable': observable, 'invariant': item})

                        #
                        if gid in Z4gid and tmp[0][-1] == ind_str:
                            body.append(tr)

                        elif tmp[0] == ind_str:
                            body.append(tr)
                        
                        count += 1
                    
                    assert len(body) > 0, 'Wrong input indicator! No such indicator in this space group!'                        

                    return {'head': head, 'body': body}
        
        else:
            raise Exception('wrong gid!')

# test:
if __name__ == '__main__':
    table = read_indicator_table(141, [0,0,0,2])
    print(table['head'])
    for tr in table['body']:
        print(tr)

# there are 2 special situations in the indicator table, some indicators has a bar overhead, some are colored blue and in a math font, e.g,
# '\bar{4}', '\textcolor{blue}{$\mathds{1}$}'
# which need to be transformed to the web language.
                    

def tex2html(string):
    if string[0] == '$':
        string = '\(' + string.strip('$') + '\)'
    
    if 'textcolor{blue}' in string:
        num = string[-4]
        string = '<span style=\"color:blue;\">' + '\(' + num + '\)' + '</span>'


    return string









