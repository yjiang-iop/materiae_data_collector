from pymongo import MongoClient
import os
import numpy as np

def insert_mat(materials, id0=1):
    # Insert the entries into the database
    client = MongoClient(port=27017)
    db = client.topmat  # database
    collection = db.materials  # collection
    collection.create_index('id', unique=True)
    for mp_id,mat in list(materials.items()):
        mat["id"] = 'CMAT%08d' % id0
        mat['mp_id'] = 1
        print(mat.keys())
        collection.insert_one(mat)
        id0 += 1
        print('inserted %s into the database as %s' % (mat['mp_id'],mat["id"]))


def update_mat(materials, id0=1):
    # Insert the entries into the database
    client = MongoClient(port=27017)
    db = client.topmat  # database
    collection = db.materials  # collection
    for mp_id, mat in list(materials.items()):
        mat0 = collection.find_one({'mp_id': mp_id})
        mat0 = mat0 if mat0 else {}
        # print(mat0)
        # mat = mat.__dict__
        # mat0.update(mat)
        if 'id' not in mat0:
            mat0['id'] = 'CMAT%08d' % id0
            id0 += 1
        collection.update_one({'mp_id': mp_id}, {"$set": mat0}, upsert=True)
        print('updated %s in the database' % mp_id, mat0['id'])

def load_data(path):
    data = np.load('%s/mat_data.npy' % path, allow_pickle=True)
    return data

if __name__ == '__main__':
    id = None
    path=os.getcwd() + '/example/SnTe'
    material = load_data(path)
    print(type(material))
    data = dict(enumerate(material.flatten(), 1))
    # insert_mat(data)
    update_mat(data)
