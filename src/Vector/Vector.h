#ifndef VECTOR // Include guard : compiles the file only once
#define VECTOR

class Vector
{
public:
    explicit Vector(double x = 0, double y = 0, double z = 0);
    ~Vector();

    double norm2() const;
    double norm() const;
    void normalize();

    double operator[](int i) const; // renvoie un double constant (qu'on ne peut pas modifier) genre a[0]
    double &operator[](int i);      // renvoie la référence vers un double pas constant (pas de const pour pouvoir le modifier) genre a[0] = 2

private:
    double data[3]; // tableau statique de 3 éléments (en précision double), qu'on remplit avec les coordonnées dans le costructeur
};

// Opérateurs définis pour les vecteurs
Vector operator+(const Vector &a, const Vector &b);
Vector operator-(const Vector &a, const Vector &b);
Vector operator-(const Vector &a);
Vector operator*(const double a, const Vector &b);
Vector operator*(const Vector &a, const double b);
Vector operator/(const Vector &a, const double b);
double dot(const Vector &a, const Vector &b);
Vector cross(const Vector &a, const Vector &b);
Vector operator*(const Vector &a, const Vector &b); // Produit terme à terme

#endif