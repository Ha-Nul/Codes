#include<iostream>
#include<string>
#include<cstring>
using namespace std;

/*
void mySort(int arr[],int n)
{
    int a=0;
    int b=0;

    while(true)
    {
        a = a + 1;
        try
        {
            if(arr[a-1] > arr[a])
            {
                int c = arr[a-1];
                int d = arr[a];
                arr[a] = c;
                arr[a-1] = d;
            }
        }
        catch(...)
        {
            //None
        }   
        for(int e = 0; e < n; ++e)
        {
            try
            {
                if(arr[e] <= arr[e+1])
                {
                    b = b + 1;
                }
                else if(arr[e] > arr[e+1])
                {
                    b = 0;
                }
            }
            catch(...)
            {
                //None
            }
        }
        if (b > n)
        {
            break;
        }
        else if (a == n)
        {
            a = 0;
        }
    }
    for (int i=0; i < n; i++)
    {
        cout << arr[i] << endl;
    }
    
}

void printArray(int arr[], int size)
{
    for (int i=0; i < size; i++)
    {
        cout << arr[i] << endl;
    }
}



int main()
{
    int arr[] = {3,1,2,5,6,7,10};
    int n = sizeof(arr)/sizeof(arr[0]);
    mySort(arr,n);
    cout << "Sorted array: \n";
    printArray(arr,n);

    return 0;
}
*/

/*
void reverse(const char* str)
{
    char arr[strlen(str) + 1];
    int a = strlen(str);
    
    for(int i=0 ; i<strlen(str) ; i++)
    {
        arr[i] = str[a-1];
        a = a - 1;
    }

    cout << arr << endl;
}

int main()
{
    const char* str = "Happy Day";
    reverse(str);
    return 0;
}
*/

/*
string asdfasdf(int x, const string &y ,int num_sentences) 
{
    string result;
    
    for (size_t i = 0; i < y.length(); ++i) {
        char currentChar = y[i];
        
        if ((currentChar >= 'A' && currentChar <= 'A' + x - 1) ||
            (currentChar >= 'a' && currentChar <= 'a' + x - 1)) {
            result += currentChar;
        } else {
            result += '#';
        }
    }

    if (num_sentences != 0)
    {
        while (num_sentences != 0)
        {
            num_sentences = num_sentences-1
            asdfasdf
        }
    }
    
    return result;
}

int main()
{
    cout << asdfasdf(4,"Hello");

    return 0;
}
*/

/*
void Test(string x)
{
    for ( int i ; i<x.length() ; i++ )
    {
        x[i] = 22;
    }

    cout << x;
}

int main()
{
    string A = "String test";
    Test(A);
}
*/

/*
int main() {
 
//Print NEWS
 
const char* sentence[] = { "HELLO CLASS", "GOOD NIGHT WORLD", "BE HAPPY" };
 
cout << (char)*(sentence[1] + 5);
 
cout << (char)(*(*(sentence + 2) + (1)));
 
cout << (char)(*(*(sentence + (1)) + 11));
 
cout << (char)(*(*sentence + (9)));
 
return 0;
}
*/

char** special_filter(const char* sentences[], int filter_offset, int num_sentences)
{
    for (int i = 0; i < num_sentences; i++) {
        const char* Nuts = sentences[i];
        for (int j = 0; j< strlen(Nuts);j++ )
        {
            int Nutss = Nuts[j];
            int Nutsss = 65+filter_offset;
            int Nutssss = 97+filter_offset;

            if (Nutsss<=Nutss && Nutss<=90)
            {
                Nutss = 35;
            }
            if (Nutssss<=Nutss && Nutss<=122)
            {
                Nutss = 35;
            }
            char stuN = Nutss;
            cout << stuN;
        } 
        //cout << strlen(Nuts) << endl;
    }
    return 0; 
}

int main()
{
    const char* sentences[] = {"Hello brus", "I am ABCDEFGHIJKLMNOPQRSTUVWXYZ"};
    int num_sentences = sizeof(sentences) / sizeof(sentences[0]);

    special_filter(sentences, 6, num_sentences);

    return 0;
}