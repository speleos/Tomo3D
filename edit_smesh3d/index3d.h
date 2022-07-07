/*
 * index3d.h based on index.h
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#ifndef _TOMO_INDEX3D_H_
#define _TOMO_INDEX3D_H_

class Index2d {
public:
    Index2d(){}
    Index2d(int i1, int i2){ i_[0]=i1; i_[1]=i2; }

    int i() const { return i_[0]; }
    int j() const { return i_[1]; }
    void set(int i1, int i2) { i_[0]=i1; i_[1]=i2; }
    
private:
    int i_[2];
};

class Index3d {
public:
    Index3d(){}
    Index3d(int i1, int i2, int i3){ i_[0]=i1; i_[1]=i2; i_[2]=i3; }

    int i() const { return i_[0]; }
    int j() const { return i_[1]; }
    int k() const { return i_[2]; }
    void set(int i1, int i2, int i3) { i_[0]=i1; i_[1]=i2; i_[2]=i3; }
    
private:
    int i_[3];
};

#endif /* _TOMO_INDEX3D_H_ */
