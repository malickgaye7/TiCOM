#include "utility.h"
//#include "tirade.h"

// ....................................................... DATA TRANSFER STUFF .......................................................


unsigned char check_container_members(Tiradata* container) { // TODO
    unsigned char container_validity = 0x01; // container is valid by default.

    if ((container -> type == "heating_axi") && ((!(container -> etheta)) || (!(container -> etheta)))) {
        container_validity = 0x00;
    }
    // TODO: DO THIS FOR ALL MEMBERS RELEVANT TO SPECIFIC TYPE (find a more efficient way of doing this?)

    return container_validity;
}

unsigned char flags_to_Tiradata(char* flags) { // TODO
    int num_flags = strlen(flags);
    Tiradata** containers = malloc(sizeof(Tiradata*) * num_flags);

    // TODO
    //  - initialize Tiradata** (pointer to array of containers) in heap
    //  - if flag contains relevant character, generate container
    //  - ^ do this as many times as there is characters in the flags string
    //  - add containers to Tiradata** array
    //  - return Tiradata** ptr
}

/************************************************************
* construct_Tiradata: Constructs a heap struct & allocates  *
*   memory based the type of data we want it to contain.    *
* Parameters:                                               *
*   char type_data[] - String indicating which data we need *
*       to store in dynamic memory                          *
*   int* input_ints - SEE COMMENT BOX BELOW (important)     *
************************************************************/
Tiradata* construct_Tiradata(char* type_data, int* input_ints) { // FINISHED (?)
    /* INPUT_INT FORMATTING:

         TYPE_DATA ............ FORMAT OF INPUT_INTS

         heating_axi: input_ints[3] = { grids.etheta size, grids.ephi size, 2d heating_axi matrix size }
         heating_radial, heating_radial_ice: input_ints[2] = { grids.etheta size, heating_radial size }
         potential: input_ints[3] = { grids.eradius size, grids.ephi size, 2d potential matrix size }
         y_i: input_ints[2] = { grids.eradius size, 2d y_i matrix size }
         lovenumber: input_ints = NULL
         Imk2: input_ints = NULL
         sensparam: input_ints[2] = { grids.eradius size, # of 'for' loop iterations }
         stress_strain: input_ints[3] = { grids.etheta size, grids.ephi size, grids.eradius size, # of 'for' loop iterations }

    */

    Tiradata* container = malloc(sizeof(Tiradata)); // initialize in heap memory
    if (!container) { return NULL; } // bad pointer, no more space
    
    // Name the structure based on its type (backend only)
    container -> type = malloc(sizeof(type_data) + 1);
    strcpy(container -> type, type_data);

    // Allocate memory based on type_data
    if (type_data == "heating_axi") {
		container -> etheta = malloc(sizeof(long double) * (*input_ints));
		container -> eradius = malloc(sizeof(long double) * (*(input_ints + 1)));
		container -> heating_axi = malloc(sizeof(long double) * (*(input_ints + 2)));
    }
    if ((type_data == "heating_radial") || (type_data == "heating_radial_ice")) {
        container -> eradius = malloc(sizeof(long double) * (*input_ints));
        container -> heating_radial = malloc(sizeof(long double) * (*(input_ints + 1)));
    }
    if (type_data == "potential") {
        container -> eradius = malloc(sizeof(long double) * (*input_ints));
        container -> ephi = malloc(sizeof(long double) * (*(input_ints + 1)));
        container -> potential = malloc(sizeof(long double) * (*(input_ints + 2)));
    }
    if (type_data == "y_i") {
        container -> eradius = malloc(sizeof(long double) * (*input_ints));
        container -> y_i_real = malloc(sizeof(long double) * (*(input_ints + 1)));
        container -> y_i_imag = malloc(sizeof(long double) * (*(input_ints + 1)));
    }
    if (type_data == "lovenumber") {
        container -> Loveh2_real = malloc(sizeof(long double));
        container -> Loveh2_imag = malloc(sizeof(long double));
        container -> Lovek2_real = malloc(sizeof(long double));
        container -> Lovek2_imag = malloc(sizeof(long double));
        container -> phaseh2 = malloc(sizeof(long double));
        container -> phasek2 = malloc(sizeof(long double));
    }
    if (type_data == "Imk2") {
        container -> Lovek2_imag = malloc(sizeof(long double));
    }
    if (type_data == "sensparam") {
        container -> eradius = malloc(sizeof(long double) * (*input_ints));
        container -> Hmu = malloc(sizeof(long double) * (*(input_ints + 1)));
        container -> hHmu = malloc(sizeof(long double) * (*(input_ints + 1)));
        container -> mu_real = malloc(sizeof(long double) * (*(input_ints + 1)));
        container -> mu_imag = malloc(sizeof(long double) * (*(input_ints + 1)));
    }
    if (type_data == "stress_strain") {
        container -> etheta = malloc(sizeof(long double) * (*input_ints));
        container -> ephi = malloc(sizeof(long double) * (*(input_ints + 1)));
        container -> eradius = malloc(sizeof(long double) * (*(input_ints + 2)));
        container -> nodes_ordered = malloc(sizeof(int) * (*(input_ints + 3)));
        container -> stress_real = malloc(sizeof(long double) * (*(input_ints + 3)));
        container -> stress_imag = malloc(sizeof(long double) * (*(input_ints + 3)));
        container -> strain_real = malloc(sizeof(long double) * (*(input_ints + 3)));
        container -> strain_imag = malloc(sizeof(long double) * (*(input_ints + 3)));
    }

    if (!check_container_members(container)) { return NULL; }

    return container; // return pointer to heap struct
}


/************************************************************
* destroy_Tiradata: Takes in a container as parameter,      *
*    and destroys it along with its relevant members.       *
* Parameters                                                *
*    Tiradata* container - struct w/ data in heap memory    *
* Return                                                    *
*    void                                                   *
************************************************************/
void destroy_Tiradata(Tiradata* container) {
    if (container -> type == "heating_axi") {
        free(container -> etheta);
        free(container -> eradius);
        free(container -> heating_axi);
    }
    if ((container -> type == "heating_radial") || (container -> type == "heating_radial_ice")) {
        free(container -> eradius);
        free(container -> heating_radial);
    }
    if (container -> type == "potential") {
        free(container -> potential);
        free(container -> ephi);
        free(container -> eradius);
    }
    if (container -> type == "y_i") {
        free(container -> eradius);
        free(container -> y_i_real);
        free(container -> y_i_imag);
    }
    if (container -> type == "lovenumber") {
        free(container -> Loveh2_imag);
        free(container -> Loveh2_real);
        free(container -> Lovek2_imag);
        free(container -> Lovek2_real);
        free(container -> phaseh2);
        free(container -> phasek2);
    }
    if (container -> type == "Imk2") {
        free(container -> Lovek2_imag);
    }
    if (container -> type == "sensparam") {
        free(container -> eradius);
        free(container -> Hmu);
        free(container -> hHmu);
        free(container -> mu_real);
        free(container -> mu_imag);
    }
    if (container -> type == "stress_strain") {
        free(container -> eradius);
        free(container -> etheta);
        free(container -> ephi);
        free(container -> stress_imag);
        free(container -> stress_real);
        free(container -> strain_imag);
        free(container -> strain_real);
        free(container -> nodes_ordered);
    }

    free(container -> type); // free name
    free(container); // free the structure itself
}

// ....................................................... UTILITY FUNCTIONS (modularized) .......................................................

// Relevant Variable Declarations
struct unitFlags uFlags;
struct matprop *layers;

long double getAvg(long double* data_ptr, int column, int numRows, int numCols) {
    int row;
    long double sum = 0;
    for (row = 0; row < numRows; row++) { sum += *(data_ptr + row + column * numRows); }
    return sum / numRows;
}

long double linearApprox(long double time, int column, long double* data_ptr, int numRows, int numCols) {
    long double t1, t2, var1, var2;
    int row; // counter
    unsigned char dataFound = OFFBYTE;

    // Find the time & variable data closest to the requested time
    for (row = 0; row < numRows - 1; row++) {
        t1 = *(data_ptr + row);
        t2 = *(data_ptr + row + 1);

        if (t1 == time) { return *(data_ptr + row + column * numRows); }
        else if (t2 == time) { return *(data_ptr + row + 1 + column * numRows); }
        else if ((t1 < time) && (t2 > time)) {
            var1 = *(data_ptr + row + column * numRows);
            var2 = *(data_ptr + row + 1 + column * numRows);
            dataFound = ONBYTE; // set to a character that can be read as a boolean true
            break;
        }
    }
    // If data not found, return average of layer data over time
    if (!dataFound && (time > t2)) { // if time out of reach
        var1 = *(data_ptr + numRows - 2 + column * numRows);
        var2 = *(data_ptr + numRows - 1 + column * numRows);
    }
    else if (!dataFound) { return getAvg(data_ptr, column, numRows, numCols); } // for edge cases
    // Return the linear estimation (used point-slope formula)
    return ((var2 - var1) / (t2 - t1)) * (time - t2) + var2;
}

/************************************************************
* layerData -- Wrapper for accessing matprop data, for both *
*       time independent & time dependent data. Give a      *
*       layer, a readable flag, along with the bit flag,    *
*       along with other necessary data (or give 0/void if  *
*       inapplicable).                                      *
* Parameters                                                *
*  int layer - layer index (starts at layer 0 -> second     *
*    column in a CSV with time dependent data)              *
*  char flag - 'r', 'd', 's', 'b', 'v'; indicates either    *
*    radius, density, shear mod, bulk mod, or visc data     *
*  char* path_flags_found - the character pointer with bits *
*    set according to a given variable is time dependent    *
*  long double time - time for which we want data; enter    *
*    0.0 if variable time independent                       *
*  int numRows - number of rows the user requested to be    *
*    read; enter 0 if variable time independent             *
*  int numCols - number of max layers + 1; enter 0 if       *
*    variable time independent                              *
*  long double* data_ptr - pointer to data in dynamic       *
*    memory; use void long double pointer if variable       *
*    is time independent                                    *
* Return (type - long double; conditional)                  *
*   1) dereferenced pointer to matprop variable if the      *
*      variable is time independent                         *
*  2a) linearly fitted point based on data corresponding to *
*      the two times closest to the requested time          *
*  2b) average of layer data for user-requested # of rows   *
************************************************************/
long double layerData(int layer, char flag, char* path_flags_found, long double time, int numRows,
                                                                int numCols, long double* data_ptr) {
    char flagList[] = "rdsbvulAg";
    long double* matprop_ptrs[9] = { &layers[layer].rad,
                                     &layers[layer].dens,
                                     &layers[layer].shear,
                                     &layers[layer].bulk,
                                     &layers[layer].visc,
                                     &layers[layer].rigidity,
                                     &layers[layer].lambda,
                                     &layers[layer].A,
                                     &layers[layer].g };
                                     // array can accommodate more/all matprop variables if necessary
    // TODO: long double & long double complex need to be placed in different arrays... rigidity & lambda
    int i;
    for (i = 0; i < 9; i++) {
        if ((i > 4 || !(*path_flags_found & (1 << i))) && (flag == flagList[i])) { // if flag & time-independent
            return *matprop_ptrs[i];
        }
        else if (flag == flagList[i]) { // if flag & time dependent
            return linearApprox(time, layer + 1, data_ptr, numRows, numCols);
        }
    }
}

/************************************************************
 * flagSearch -- Reads line from input file and returns an  *
 *      unsigned character that indicates if any pathname   *
 *      flags have been found in that line.                 *
 *                                                          *
 * Parameters                                               *
 *  line -- any given line from the input file              *
 ************************************************************/
// FINISHED
unsigned char flagSearch(char line[]) {
    unsigned char path_flags = OFFBYTE; // bit representation initialized to 00000000
    int i; // counter

    char flags[5][4] = { "rad", "den", "she", "bul", "vis" }; // flags to find
    for (i = 0; i < 5; i++) {
        if (strstr(line, flags[i]) != NULL) { path_flags |= 1 << i; }
    }
    /* If flag found, i-th bit from right is set to 1
            "rad" flag found -> xxxxxxx1
            "den" flag found -> xxxxxx1x
            "she" flag found -> xxxxx1xx
            "bul" flag found -> xxxx1xxx
            "vis" flag found -> xxx1xxxx
        ... if multiple flags found, then multiple bits are set to 1
    */
    return path_flags;
}

void stitchLine(char line[], char* path_flags_found) {
    int i;
    for (i = 0; i < 5; i++) {
        if (i != 1 && !(*path_flags_found & (1 << i))) { strcat(line, "%Le "); }
        else if (!(*path_flags_found & (1 << i))) { strcat(line, "%Lf "); }
    }
    line[strlen(line)-1] = '\0'; // change last whitespace character to a string terminator
}

void setUnits(char* label, char* path_flags_found) {
    unsigned char* flag_ptrs[5] = { &uFlags.tUnits,
                                    &uFlags.radUnits,
                                    &uFlags.densUnits,
                                    &uFlags.shearbulkUnits,
                                    &uFlags.viscUnits };
    int i;
    for (i = 0; i < 5; i++) { if (*flag_ptrs[i] == '\0') { *flag_ptrs[i] = OFFBYTE; } } // set flags to empty if null
    char unitStrings[25][12] = { "(years)", "(months)", "(weeks)", "(days)", "(hours)", "(min)", "(s)", // Time flags (0 - 6)
                                 "(km)", "(m)", "(cm)", "(mi)", "(yd)", "(ft)", "(in)", // Radius flags (7 - 13)
                                 "(kg/m^3)", "(g/cm^3)", "(lb/ft^3)", // Density flags (14 - 16)
                                 "(GPa)", "(Pa)", "(ksi)", "(psi)", // Shear & Bulk modulus flags (17 - 20)
                                 "(Ns/m^2)", "(Pa s)", "(kg/ms)", "(lb s/ft^2)" }; // Viscosity flags (21 - 24; 21-23 are equivalent)

    for (i = 0; i < 24; i++) {
        if (strstr(label, unitStrings[i]) != NULL) {
            switch (i) {
                case 0 ... 6: // if time unit flag
                    *flag_ptrs[0] |= 1 << i;
                    printf("\tTime units detected in file: %s\n", unitStrings[i]);
                    break;
                case 7 ... 13: // if radius unit flag
                    *flag_ptrs[1] |= 1 << (i - 7);
                    printf("\tRadius units detected in file: %s\n", unitStrings[i]);
                    break;
                case 14 ... 16: // if density unit flag
                    *flag_ptrs[2] |= 1 << (i - 14);
                    printf("\tDensity units detected in file: %s\n", unitStrings[i]);
                    break;
                case 17 ... 20: // if shear/bulk unit flag
                    *flag_ptrs[3] |= 1 << (i - 17);
                    printf("\tShear/bulk modulus units detected in file: %s\n", unitStrings[i]);
                    break;
                case 21 ... 23: // if viscosity unit flag (unit 1)
                    *flag_ptrs[4] |= 1 << 0;
                    printf("\tViscosity units detected in file: %s\n", unitStrings[i]);
                    break;
                default: // if viscosity unit flag (unit 2)
                    *flag_ptrs[4] |= 1 << 1;
                    printf("\tViscosity units detected in file: %s\n", unitStrings[24]);
                    break;
            }
        }
    }
}

void convertUnits(long double* data_ptr, char* path_flags_found, int numRows, int numCols) {
    int i, j, k;

    if (!(uFlags.tUnits & (1 << 0))) { // if unit is not in years
        long double tConversions[6] = { 1.0/12.0, 1.000/52.143, 1.0/365.0, 1.0/8670.0, 1.0/525600.0, 1.0/31536000.0 };
                                    // mon to yr - wk to yr - day to yr - hr to yr - min to yr - s to yr
        for (i = 1; i < 7; i++) { // check each non-default bit
            if ((uFlags.tUnits & (1 << i))) { // if bit is set (unit found)
                printf("\tConverting time to years...\n");
                for (j = 0; j < numRows; j++) {
                    *(data_ptr + j) *= tConversions[i - 1];
                }
                break;
            }
        }
        uFlags.tUnits = OFFBYTE;
        uFlags.tUnits |= (1 << 0); // set flag to default
    }
    if ((*path_flags_found & (1 << 0)) && !(uFlags.radUnits & (1 << 0))) { // if radius unit is not in km
        long double radConversions[6] = { 1.0/1000.0, 1.0/100000.0, 1.0/1.609, 1.0/1094.0, 1.0/3281.0, 1.0/39370.0 };
                                        // m to km - cm to km - mi to km - yrd to km - ft to km - in to km
        for (i = 1; i < 7; i++) { // check each non-default bit
            if ((uFlags.radUnits & (1 << i))) { // if bit is set (unit found)
                printf("\tConverting radius to kilometers...\n");
                for (k = 1; k < numCols; k++) {
                    for (j = 0; j < numRows; j++) {
                        *(data_ptr + numRows * k + j) *= radConversions[i - 1];
                    }
                }
                break;
            }
        }
        uFlags.radUnits = OFFBYTE;
        uFlags.radUnits |= (1 << 0); // set flag to default
    }
    if ((*path_flags_found & (1 << 1)) && !(uFlags.densUnits & (1 << 0))) { // if density unit is not in kg/m^3
        long double densConversions[2] = { 1000.0, 16.0185 };
                                        // g/cm^3 to kg/m^3 - lb/ft^3 to kg/m^3
        for (i = 1; i < 3; i++) { // check each non-default bit
            if ((uFlags.densUnits & (1 << i))) { // if bit is set (unit found)
                printf("\tConverting density to kilograms per cubic meter...\n");
                for (k = 1; k < numCols; k++) {
                    for (j = 0; j < numRows; j++) {
                        *(data_ptr + numRows * k + j) *= densConversions[i - 1];
                    }
                }
                break;
            }
        }
        uFlags.densUnits = OFFBYTE;
        uFlags.densUnits |= (1 << 0); // set flag to default
    }
    if ((*path_flags_found & (1 << 2) | *path_flags_found & (1 << 3)) && !(uFlags.shearbulkUnits & (1 << 1))) {
        // if shear/bulk mod. unit is not in Pa
        long double modConversions[3] = { 1000.0, 6894757.2931783, 6894.76 };
                                        // GPa to Pa - ksi to Pa - psi to Pa
        int index;
        if (uFlags.shearbulkUnits & (1 << 0)) { index = 0; }
        if (uFlags.shearbulkUnits & (1 << 2)) { index = 1; }
        if (uFlags.shearbulkUnits & (1 << 3)) { index = 2; }

        printf("\tConverting bulk & shear modulus to Pascals...\n");
        for (k = 1; k < numCols; k++) {
            for (j = 0; j < numRows; j++) {
                *(data_ptr + numRows * k + j) *= modConversions[index - 1];
            }
        }

        uFlags.shearbulkUnits = OFFBYTE;
        uFlags.shearbulkUnits |= (1 << 1); // set flag to default
    }
    if ((*path_flags_found & (1 << 4)) && !(uFlags.viscUnits & (1 << 0))) { // if density unit is not in kg/m^3
        printf("\tConverting viscosity to Pascal-seconds...\n");
        for (k = 1; k < numCols; k++) {
            for (j = 0; j < numRows; j++) { *(data_ptr + numRows * k + j) *= 47.880259; } // lb s/ft^2 to Pa*s conversion
        }
        uFlags.viscUnits = OFFBYTE;
        uFlags.viscUnits |= (1 << 0); // set flag to default
    }
}

// 
long double* readCSVData(int numRows, int numLayers, FILE* fp, char* flag, char* path_flags_found) {
    int i = 0;
    int j = 0;
    int numCols = numLayers + 1;

    char row[80]; // row or line from input file with time-independent data
    long double *tdepData = (long double *)malloc(numCols * numRows * sizeof(long double));

    // Read column headers (flagging only; no heap memory studio)
    fgets(row, sizeof(row), fp); // get labels
    char* labels = strtok(row, ",;\t\n");
    while (labels != NULL) {
        setUnits(labels, path_flags_found);
        labels = strtok(NULL, ",;\t\n"); // next token
    }
    
    // Set default units if needed
    if (uFlags.tUnits == OFFBYTE) {
        printf("\tTime unit unclear; defaulting to years\n");
        uFlags.tUnits |= (1 << 0);
        // TODO: change default to seconds (SI units)
    }
    if ((uFlags.radUnits == OFFBYTE) && (*path_flags_found & (1 << 0))) {
        printf("\tRadius unit unclear; defaulting to kilometers\n");
        uFlags.radUnits |= (1 << 0);
    }
    if ((uFlags.densUnits == OFFBYTE) && (*path_flags_found & (1 << 1))) {
        printf("\tDensity unit unclear; defaulting to kg/m^3\n");
        uFlags.densUnits |= (1 << 0);
    }
    if ((uFlags.shearbulkUnits == OFFBYTE) && ((*path_flags_found & (1 << 2) | *path_flags_found & (1 << 3)))) {
        printf("\tShear/bulk modulus unit unclear; defaulting to Pa\n");
        uFlags.shearbulkUnits |= (1 << 1);
    }
    if ((uFlags.viscUnits == OFFBYTE) && (*path_flags_found & (1 << 4))) {
        printf("\tViscosity unit unclear; defaulting to Pa*s\n");
        uFlags.viscUnits |= (1 << 0);
    }

    // Read & store numerical data in dynamic memory
    while (fgets(row, sizeof(row), fp) && (j < numRows)) {
        i = 0; // Go to the leftmost column
        char* dataPoint = strtok(row, ",;\t\n");

        while (dataPoint != NULL) {
            char formatter[25] = "";
            stitchLine(formatter, flag);
            sscanf(dataPoint, formatter, (tdepData + numRows * i + j));
            dataPoint = strtok(NULL, ",;\t\n");
            i++; // Move to the next column (scan right)
        }
        j++; // Move to the next row (scan down)
    }
    convertUnits(tdepData, path_flags_found, numRows, numCols);
    return tdepData; // Return pointer to dynamic memory
}