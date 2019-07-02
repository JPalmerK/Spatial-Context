function [ang1 ang2 ang3] = sssTriangle(a, b, c)

ang1=acosd((a^2 + b^2 - c^2)/(2*a*b));
ang2=acosd((b^2 + c^2 - a^2)/(2*b*c));
ang3=acosd((c^2 + a^2 - b^2)/(2*a*c));

end