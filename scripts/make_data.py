from src import angles as ang
from math import ceil
import matplotlib.pyplot as plt
import json

def possDend(per,max_pre_rng=2):
    candidates = [[2*k+1,(2**per-1)*2**(n)] for n in range(1,max_pre_rng) for k in range(0,ceil(((2**per-1)*2**(n-1)-1)/2))]
    working={}
    
    # X = []
    # Y = []
    for key in candidates: 
        theta=ang.Angle(key[0],key[1])
        theta.itin_to_rat()
        angle = f"{key[0]}/{key[1]}"
        theta.itin_to_rat('^')
        working[angle]={'ks':theta.ks, 'ks_per':theta.per_len, 'start_index_per':theta.start_index_per, 'itin':theta.itin, 'itin_per':theta.itin_per_len, 'rat_func':theta.rat_func}
    return working
    #     theta.assoc_lambda()
    #     l=theta.lam
    #     try:
    #         l.evalf()
    #     except:
    #         # print(f"theta {l} cannot be evaluated")
    #         continue
    #     try:
    #         if ( (complex(l).imag!=0)  ):
    #             # working.append(key)
    #             working[angle]={'ks':theta.ks, 'ks_per':theta.per_len, 'start_index_per':theta.start_index_per, 'itin':theta.itin, 'itin_per':theta.itin_per_len, 'rat_func':theta.rat_func}

    #             X.append(complex(l).real)
    #             if (complex(l).imag<0):
    #                 Y.append(-complex(l).imag)
    #                 working[angle].update({'lambda': {"real":complex(l).real,"imag":-complex(l).imag}})
    #             else:
    #                 Y.append(complex(l).imag)
    #                 working[angle].update({'lambda':  {"real":complex(l).real,"imag":complex(l).imag}})
    #     except:
    #         # print(f"theta {l} does not work")
    #         continue
    # return working,X,Y


if __name__ == "__main__":
    # print("Start...\n")
    # perThree = possDend(3,11)
    # # perThree,x_three_coord,y_three_coord = possDend(3,2)
    # print("...finished period 3...\n")
    # perFive = possDend(5,7)
    # # perFive,x_five_coord,y_five_coord = possDend(5,2)
    # print("...finished period 5...\n")

    # perEleven = possDend(11,2)
    # print("...finished period 11...\n")
    # angle_data={**perThree,**perFive,**perEleven}
    # print("...merged the data...\n")

    # with open("data/data.json","w") as file_data:
    #     json.dump(angle_data, file_data)
    
    # print("Done.")

    
    # sorted_per={"1":{},"3":{},"5":{},"11":{}}

    # with open("data/data.json","r") as stuff:
    #     data_json=json.load(stuff)
    #     for theta in data_json.keys():
    #         period_len=f"{data_json[theta]['ks_per']}"
    #         sorted_per[period_len].update({theta: data_json[theta]})
    #     sorted_data = open("data/sorted_data.json","w")
    #     json.dump(sorted_per,sorted_data)