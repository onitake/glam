unsigned int sky_d_width = 32;
unsigned int sky_d_height = 32;
unsigned char sky_d_samples[] = { 0x50,0x67,0x83,0x50,0x66,0x83,0x4f,0x66,0x83,0x4e,0x65,0x82,0x4c,0x63,0x81,0x47,0x5f,0x7e,0x44,0x5c,0x7b,0x40,0x58,0x77,0x3b,0x52,0x72,0x39,0x51,0x71,0x37,0x4e,0x6e,0x35,0x4c,0x6c,0x32,0x49,0x6a,0x30,0x47,0x67,0x2f,0x46,0x66,0x2e,0x44,0x64,0x2e,0x43,0x63,0x2d,0x42,0x62,0x2c,0x41,0x60,0x2c,0x40,0x5f,0x2b,0x3f,0x5f,0x2b,0x3e,0x5d,0x2b,0x3e,0x5d,0x2c,0x3e,0x5c,0x2c,0x3d,0x5b,0x2c,0x3d,0x59,0x2c,0x3d,0x59,0x2a,0x3c,0x57,0x29,0x3a,0x56,0x27,0x39,0x54,0x27,0x38,0x52,0x26,0x37,0x51,0x58,0x6e,0x87,0x56,0x6c,0x87,0x54,0x6b,0x86,0x52,0x6a,0x86,0x50,0x68,0x85,0x4b,0x63,0x81,0x47,0x5f,0x7e,0x42,0x5a,0x7a,0x3d,0x56,0x75,0x3c,0x54,0x74,0x39,0x51,0x72,0x37,0x4e,0x6f,0x35,0x4b,0x6b,0x32,0x49,0x6a,0x32,0x48,0x68,0x31,0x47,0x67,0x30,0x46,0x66,0x30,0x45,0x64,0x2e,0x43,0x63,0x2d,0x42,0x62,0x2d,0x41,0x61,0x2c,0x40,0x5f,0x2c,0x40,0x5f,0x2d,0x3f,0x5e,0x2c,0x3e,0x5d,0x2d,0x3e,0x5b,0x2f,0x40,0x5c,0x2d,0x3e,0x5b,0x2c,0x3d,0x58,0x29,0x3b,0x56,0x28,0x3a,0x55,0x28,0x39,0x53,0x5b,0x71,0x8c,0x5b,0x70,0x8b,0x5a,0x6f,0x8b,0x58,0x6f,0x8b,0x55,0x6c,0x89,0x51,0x68,0x85,0x4d,0x64,0x82,0x46,0x5d,0x7d,0x48,0x5f,0x7d,0x48,0x5d,0x7b,0x40,0x56,0x76,0x38,0x50,0x71,0x36,0x4d,0x6e,0x35,0x4c,0x6d,0x33,0x4a,0x6b,0x32,0x48,0x69,0x32,0x48,0x68,0x31,0x47,0x67,0x31,0x46,0x66,0x2f,0x44,0x64,0x2e,0x43,0x63,0x2e,0x42,0x62,0x2e,0x42,0x61,0x2e,0x40,0x60,0x2d,0x40,0x5f,0x2f,0x40,0x5e,0x30,0x42,0x5e,0x3b,0x4b,0x64,0x2d,0x3e,0x5a,0x2b,0x3c,0x57,0x2a,0x3c,0x57,0x29,0x3a,0x56,0x62,0x77,0x91,0x60,0x76,0x90,0x5e,0x75,0x90,0x5d,0x75,0x90,0x5a,0x72,0x8e,0x57,0x6e,0x8b,0x52,0x68,0x86,0x4a,0x62,0x81,0x53,0x67,0x82,0x52,0x65,0x80,0x45,0x5c,0x7b,0x3b,0x52,0x73,0x38,0x4f,0x70,0x37,0x4e,0x6e,0x34,0x4b,0x6c,0x34,0x4a,0x6b,0x34,0x4a,0x69,0x33,0x49,0x69,0x32,0x48,0x68,0x31,0x46,0x66,0x31,0x45,0x65,0x32,0x47,0x65,0x30,0x45,0x64,0x30,0x44,0x62,0x30,0x43,0x62,0x30,0x43,0x61,0x32,0x44,0x60,0x32,0x45,0x61,0x2f,0x42,0x5e,0x2d,0x40,0x5b,0x2c,0x3d,0x59,0x2a,0x3c,0x57,0x63,0x79,0x94,0x62,0x78,0x95,0x62,0x79,0x95,0x63,0x79,0x95,0x60,0x77,0x92,0x5c,0x72,0x8e,0x59,0x6f,0x8b,0x53,0x69,0x85,0x59,0x6d,0x87,0x4d,0x63,0x81,0x57,0x69,0x82,0x4e,0x61,0x7d,0x3b,0x52,0x72,0x37,0x4e,0x70,0x35,0x4c,0x6e,0x35,0x4c,0x6d,0x35,0x4c,0x6d,0x35,0x4b,0x6b,0x33,0x49,0x6a,0x33,0x49,0x69,0x32,0x48,0x67,0x33,0x48,0x67,0x33,0x47,0x66,0x32,0x46,0x65,0x34,0x47,0x66,0x33,0x47,0x64,0x33,0x46,0x62,0x34,0x46,0x62,0x31,0x44,0x60,0x2f,0x42,0x5e,0x2e,0x3f,0x5b,0x2b,0x3c,0x58,0x67,0x7e,0x9a,0x67,0x7e,0x9a,0x67,0x7e,0x9b,0x6a,0x7f,0x9b,0x62,0x79,0x95,0x5c,0x72,0x90,0xa4,0xac,0xba,0x91,0x9c,0xab,0x7f,0x8b,0x9e,0x7b,0x86,0x99,0xab,0xb0,0xb9,0x60,0x70,0x87,0x3e,0x54,0x75,0x3e,0x55,0x75,0x38,0x4f,0x71,0x37,0x4f,0x70,0x37,0x4e,0x6e,0x36,0x4d,0x6d,0x35,0x4b,0x6b,0x34,0x4b,0x6b,0x36,0x4b,0x6b,0x34,0x49,0x68,0x34,0x49,0x68,0x35,0x4a,0x68,0x3e,0x50,0x6b,0x38,0x4b,0x68,0x38,0x4c,0x67,0x51,0x5e,0x73,0x34,0x47,0x62,0x31,0x45,0x60,0x2f,0x41,0x5c,0x2a,0x3c,0x59,0x6d,0x84,0xa1,0x6d,0x85,0xa1,0x6f,0x88,0xa3,0x71,0x87,0xa3,0xa3,0xae,0xbd,0xb4,0xbb,0xc6,0xb5,0xbc,0xc5,0x9e,0xa5,0xb2,0xb6,0xba,0xc2,0xc1,0xc4,0xca,0xca,0xcc,0xd0,0xa3,0xa9,0xb4,0x85,0x8f,0x9f,0x4f,0x62,0x7e,0x3b,0x52,0x73,0x3a,0x51,0x72,0x39,0x50,0x71,0x38,0x4e,0x6f,0x37,0x4e,0x6f,0x37,0x4c,0x6c,0x36,0x4b,0x6c,0x37,0x4c,0x6c,0x39,0x4d,0x6b,0x3e,0x51,0x6c,0x3a,0x4d,0x6a,0x3b,0x4d,0x6a,0x3e,0x50,0x6b,0x5a,0x67,0x7a,0x4e,0x5b,0x72,0x34,0x47,0x62,0x2d,0x40,0x5e,0x2c,0x3e,0x5b,0x72,0x8b,0xa7,0x73,0x8c,0xa8,0x75,0x8d,0xa9,0x87,0x99,0xb1,0xa6,0xb1,0xbf,0xd1,0xd3,0xd8,0xd6,0xd8,0xdc,0xcd,0xd0,0xd5,0xc4,0xc6,0xcc,0xbc,0xc1,0xc9,0xbf,0xc3,0xca,0xcd,0xd0,0xd4,0xba,0xbe,0xc5,0x61,0x71,0x89,0x3f,0x56,0x78,0x58,0x69,0x83,0x40,0x57,0x76,0x3e,0x54,0x74,0x3d,0x52,0x72,0x3a,0x50,0x70,0x38,0x4d,0x6d,0x37,0x4d,0x6d,0x3a,0x4e,0x6c,0x3a,0x4f,0x6d,0x3b,0x4f,0x6d,0x42,0x55,0x6f,0x4a,0x5b,0x74,0x6e,0x77,0x86,0x7d,0x82,0x8e,0x32,0x45,0x63,0x2e,0x42,0x5f,0x2d,0x40,0x5d,0x78,0x91,0xad,0x82,0x99,0xb2,0x7c,0x95,0xaf,0x83,0x98,0xb1,0x93,0xa2,0xb5,0xd3,0xd6,0xda,0xe0,0xe0,0xe3,0xde,0xdf,0xe1,0xdc,0xdd,0xe0,0xcd,0xd0,0xd5,0xbf,0xc4,0xcc,0xb6,0xbc,0xc4,0xb2,0xb8,0xc1,0x49,0x5f,0x7d,0x45,0x5c,0x7c,0x54,0x67,0x83,0x48,0x5d,0x7b,0x3d,0x54,0x74,0x3d,0x54,0x74,0x3a,0x50,0x70,0x3a,0x50,0x70,0x3c,0x51,0x70,0x39,0x4d,0x6e,0x3c,0x51,0x6f,0x4d,0x5d,0x77,0x45,0x57,0x71,0x46,0x57,0x71,0x8d,0x94,0x9e,0x6f,0x77,0x87,0x34,0x49,0x66,0x31,0x44,0x62,0x2e,0x42,0x5f,0x84,0x9c,0xb5,0x8c,0xa2,0xba,0xa2,0xb1,0xc4,0xc0,0xc8,0xd3,0xb1,0xbc,0xc9,0xca,0xce,0xd6,0xdd,0xdf,0xe3,0xe6,0xe6,0xe8,0xda,0xdc,0xdf,0xc8,0xcd,0xd3,0xb0,0xb6,0xc0,0xa3,0xaa,0xb5,0x9a,0xa3,0xb0,0x81,0x8e,0x9e,0x43,0x5c,0x7d,0x46,0x5d,0x7d,0x59,0x6b,0x84,0x74,0x80,0x94,0x49,0x5d,0x7b,0x3a,0x51,0x72,0x3b,0x51,0x70,0x3c,0x52,0x71,0x3b,0x50,0x6f,0x38,0x4e,0x6d,0x3a,0x4f,0x6d,0x3c,0x51,0x6e,0x48,0x5a,0x74,0x97,0x9a,0xa4,0x5a,0x66,0x7b,0x36,0x4a,0x67,0x31,0x45,0x64,0x30,0x42,0x60,0x8a,0xa3,0xbb,0x8f,0xa5,0xbd,0xb2,0xbf,0xcf,0xe8,0xe8,0xeb,0xcf,0xd6,0xdd,0xb1,0xbb,0xc9,0xb3,0xbc,0xc7,0xbd,0xc4,0xce,0xb6,0xbe,0xc7,0xab,0xb3,0xbe,0xa1,0xaa,0xb5,0x9a,0xa3,0xaf,0x94,0x9e,0xac,0x89,0x94,0xa4,0x6e,0x7d,0x92,0x57,0x6b,0x88,0x54,0x68,0x84,0x6a,0x77,0x8d,0x5d,0x6d,0x86,0x3d,0x53,0x74,0x38,0x50,0x71,0x3d,0x54,0x73,0x3d,0x52,0x71,0x3b,0x50,0x70,0x3a,0x4f,0x6e,0x3b,0x50,0x6f,0x43,0x58,0x74,0x6d,0x78,0x87,0x51,0x5f,0x76,0x38,0x4b,0x68,0x32,0x46,0x65,0x30,0x43,0x62,0x8e,0xa6,0xbf,0x93,0xaa,0xc3,0xac,0xbc,0xcf,0xd2,0xd7,0xdf,0xe0,0xe4,0xe9,0xce,0xd4,0xdb,0xd1,0xd6,0xdc,0xb7,0xbf,0xca,0xa7,0xb0,0xbb,0xa2,0xaa,0xb5,0x9e,0xa7,0xb1,0x9e,0xa7,0xb2,0x99,0xa2,0xae,0x98,0xa1,0xad,0x9f,0xa6,0xb2,0xb5,0xba,0xc3,0xa5,0xab,0xb6,0xa1,0xa8,0xb3,0x6b,0x79,0x8d,0x49,0x5e,0x7d,0x40,0x56,0x75,0x42,0x58,0x77,0x45,0x5a,0x78,0x3f,0x54,0x74,0x3f,0x53,0x72,0x3c,0x52,0x70,0x40,0x54,0x72,0x64,0x70,0x83,0x53,0x62,0x79,0x39,0x4d,0x6a,0x34,0x47,0x66,0x30,0x44,0x63,0x93,0xac,0xc6,0x96,0xae,0xc8,0xc2,0xcd,0xdb,0xf7,0xf7,0xf8,0xf5,0xf6,0xf7,0xe0,0xe4,0xe7,0xca,0xd1,0xd8,0xb8,0xc1,0xcc,0xae,0xb7,0xc2,0xa9,0xb2,0xbc,0xa4,0xac,0xb8,0xa3,0xab,0xb6,0xa1,0xa9,0xb4,0xac,0xb3,0xbd,0xb5,0xbb,0xc2,0xcf,0xd2,0xd6,0xbe,0xc3,0xca,0xba,0xbe,0xc7,0x83,0x8d,0x9f,0x48,0x5d,0x7d,0x58,0x6a,0x84,0x5e,0x6e,0x86,0x44,0x5a,0x79,0x44,0x59,0x77,0x3e,0x54,0x73,0x3d,0x52,0x71,0x3c,0x51,0x6f,0x6b,0x77,0x88,0x62,0x6e,0x82,0x3a,0x4e,0x6c,0x35,0x49,0x67,0x31,0x45,0x63,0xa2,0xb9,0xd1,0xc6,0xd2,0xde,0xe4,0xe8,0xec,0xfd,0xfd,0xfd,0xff,0xff,0xff,0xf5,0xf6,0xf6,0xd2,0xd8,0xde,0xbe,0xc5,0xcf,0xb7,0xc0,0xca,0xb1,0xba,0xc4,0xac,0xb4,0xbe,0xab,0xb3,0xbe,0xb5,0xbc,0xc5,0xcb,0xcf,0xd3,0xd2,0xd4,0xd8,0xd2,0xd5,0xd9,0xc1,0xc6,0xcc,0xbd,0xc2,0xc9,0x8d,0x97,0xa7,0x58,0x6b,0x85,0x62,0x74,0x8d,0x66,0x75,0x8c,0x4c,0x60,0x7d,0x46,0x5b,0x7a,0x41,0x57,0x76,0x3d,0x53,0x72,0x3a,0x4f,0x70,0x47,0x5b,0x77,0x5b,0x68,0x7e,0x3f,0x51,0x6e,0x36,0x4a,0x69,0x32,0x46,0x65,0xbe,0xce,0xdf,0xe3,0xe8,0xed,0xd2,0xdb,0xe4,0xf1,0xf2,0xf3,0xfe,0xfe,0xfe,0xf3,0xf4,0xf5,0xdc,0xe1,0xe5,0xc7,0xce,0xd6,0xc2,0xca,0xd1,0xba,0xc1,0xcb,0xb2,0xba,0xc5,0xb4,0xbc,0xc6,0xb3,0xba,0xc4,0xc3,0xc9,0xcf,0xdd,0xde,0xe0,0xcf,0xd2,0xd7,0xbc,0xc1,0xc9,0xb8,0xbc,0xc5,0x8c,0x96,0xa6,0x58,0x6c,0x87,0x52,0x67,0x83,0x4e,0x63,0x80,0x47,0x5d,0x7d,0x43,0x5a,0x7a,0x41,0x57,0x77,0x3d,0x53,0x73,0x3a,0x51,0x70,0x38,0x4e,0x6d,0x3b,0x50,0x6e,0x3a,0x4f,0x6e,0x37,0x4b,0x69,0x3a,0x4c,0x69,0xb2,0xc8,0xdf,0xc5,0xd4,0xe4,0xce,0xda,0xe6,0xed,0xee,0xf1,0xfa,0xfa,0xfb,0xe4,0xe7,0xeb,0xe1,0xe5,0xe9,0xd0,0xd5,0xdb,0xc7,0xce,0xd5,0xc0,0xc7,0xd0,0xbe,0xc6,0xce,0xb9,0xc1,0xca,0xb3,0xbb,0xc5,0xb1,0xb8,0xc2,0xdb,0xdc,0xdf,0xd8,0xda,0xde,0xda,0xdb,0xde,0xa1,0xa9,0xb6,0x8a,0x95,0xa6,0x57,0x6b,0x87,0x50,0x66,0x84,0x4e,0x63,0x81,0x49,0x5f,0x7f,0x44,0x5b,0x7a,0x42,0x58,0x78,0x3d,0x53,0x74,0x39,0x50,0x70,0x37,0x4d,0x6e,0x39,0x4e,0x6d,0x3a,0x4f,0x6e,0x36,0x4b,0x69,0x43,0x54,0x6f,0xa5,0xc1,0xde,0xbd,0xd0,0xe4,0xef,0xf1,0xf4,0xf4,0xf5,0xf7,0xeb,0xee,0xf0,0xe6,0xe9,0xed,0xf4,0xf4,0xf5,0xd8,0xdd,0xe3,0xcd,0xd3,0xd9,0xcb,0xd0,0xd7,0xc9,0xcf,0xd6,0xc2,0xc8,0xd1,0xbb,0xc2,0xcb,0xb5,0xbd,0xc7,0xc0,0xc5,0xcd,0xd3,0xd5,0xda,0xdf,0xdf,0xe1,0x83,0x90,0xa4,0x6f,0x80,0x98,0x53,0x69,0x87,0x4d,0x64,0x83,0x4f,0x65,0x82,0x49,0x60,0x7f,0x41,0x59,0x7a,0x3f,0x56,0x77,0x3d,0x54,0x74,0x39,0x50,0x71,0x37,0x4d,0x6d,0x38,0x4d,0x6d,0x39,0x4f,0x6d,0x39,0x4d,0x6c,0x3b,0x4e,0x6c,0xaa,0xc5,0xe2,0xe5,0xeb,0xf1,0xfd,0xfd,0xfe,0xf9,0xfa,0xfa,0xe0,0xe5,0xeb,0xe8,0xea,0xee,0xf5,0xf6,0xf5,0xe6,0xe9,0xec,0xd6,0xda,0xdf,0xca,0xd0,0xd8,0xcf,0xd3,0xd9,0xc8,0xce,0xd5,0xc9,0xcf,0xd5,0xc8,0xce,0xd5,0xb2,0xba,0xc5,0xb2,0xba,0xc4,0xc8,0xca,0xd1,0x80,0x90,0xa6,0x7b,0x89,0x9f,0x56,0x6d,0x8a,0x50,0x67,0x86,0x4c,0x64,0x82,0x4c,0x62,0x81,0x47,0x5d,0x7d,0x3f,0x55,0x78,0x3b,0x53,0x73,0x3a,0x51,0x71,0x38,0x4e,0x6f,0x3a,0x4f,0x6f,0x39,0x4e,0x6e,0x38,0x4d,0x6c,0x39,0x4d,0x6b,0xbc,0xd1,0xe7,0xf9,0xfa,0xfb,0xff,0xff,0xff,0xff,0xff,0xff,0xf3,0xf4,0xf6,0xdd,0xe2,0xe9,0xec,0xef,0xf1,0xef,0xf1,0xf3,0xdb,0xe0,0xe4,0xcc,0xd2,0xda,0xc8,0xcd,0xd6,0xcb,0xd0,0xd7,0xcc,0xd1,0xd7,0xca,0xcf,0xd6,0xbd,0xc4,0xcd,0xc4,0xca,0xd1,0xb8,0xc0,0xca,0xbe,0xc5,0xcc,0xd1,0xd2,0xd7,0xc8,0xcc,0xd2,0xaa,0xb1,0xbd,0x75,0x85,0x9a,0xa1,0xa9,0xb4,0x7e,0x89,0x9b,0x45,0x5b,0x7b,0x3d,0x55,0x76,0x3a,0x52,0x72,0x3b,0x52,0x72,0x3c,0x52,0x71,0x3d,0x52,0x71,0x3a,0x4f,0x6e,0x3a,0x4e,0x6d,0xd3,0xe0,0xef,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xef,0xf2,0xf5,0xe0,0xe6,0xeb,0xe3,0xe7,0xec,0xe2,0xe6,0xeb,0xe4,0xe7,0xeb,0xd4,0xd9,0xde,0xc6,0xcd,0xd5,0xc4,0xcb,0xd3,0xc0,0xc8,0xd1,0xb6,0xbe,0xc9,0xb6,0xbe,0xc8,0xc9,0xce,0xd6,0xcc,0xd1,0xd6,0xca,0xcf,0xd5,0xd1,0xd3,0xd9,0xe5,0xe5,0xe7,0xd8,0xd9,0xdd,0xb1,0xb7,0xc1,0xbb,0xbe,0xc5,0xa9,0xae,0xb7,0x4b,0x61,0x7f,0x3e,0x55,0x77,0x3d,0x54,0x74,0x3e,0x55,0x74,0x3e,0x54,0x73,0x3f,0x54,0x72,0x3b,0x50,0x6f,0x5b,0x6a,0x7f,0xd9,0xe5,0xf3,0xfb,0xfc,0xfd,0xff,0xff,0xff,0xf9,0xf9,0xfa,0xf3,0xf3,0xf6,0xe9,0xee,0xf1,0xd8,0xdf,0xe7,0xdb,0xe1,0xe8,0xe6,0xe9,0xee,0xdb,0xdf,0xe4,0xca,0xd0,0xd9,0xbd,0xc6,0xd0,0xbe,0xc6,0xd1,0xbd,0xc5,0xcf,0xb6,0xbf,0xc9,0xcb,0xd0,0xd7,0xd7,0xdb,0xde,0xda,0xdc,0xe0,0xc6,0xca,0xd1,0xdc,0xde,0xe1,0xdd,0xdf,0xe1,0xc5,0xc9,0xcf,0x9a,0xa2,0xaf,0x76,0x83,0x97,0x51,0x65,0x81,0x3e,0x57,0x78,0x3f,0x56,0x77,0x40,0x56,0x75,0x40,0x56,0x75,0x44,0x5a,0x77,0x6d,0x7b,0x8d,0xa4,0xa8,0xaf,0xf2,0xf6,0xfa,0xff,0xff,0xff,0xff,0xff,0xff,0xfe,0xfe,0xfe,0xff,0xff,0xff,0xf3,0xf5,0xf7,0xd5,0xde,0xe8,0xd4,0xdc,0xe4,0xda,0xe0,0xe6,0xd4,0xd9,0xe1,0xc6,0xce,0xd8,0xbb,0xc6,0xd1,0xbc,0xc6,0xd1,0xbd,0xc6,0xd0,0xbc,0xc5,0xcf,0xd6,0xda,0xdf,0xe1,0xe4,0xe6,0xe8,0xe9,0xea,0xde,0xdf,0xe2,0xd9,0xdb,0xde,0xe5,0xe6,0xe7,0xda,0xdb,0xdd,0xbc,0xc1,0xc7,0x6f,0x7f,0x95,0x4a,0x61,0x80,0x43,0x5b,0x7c,0x40,0x58,0x78,0x43,0x59,0x7a,0x46,0x5b,0x79,0x80,0x8a,0x99,0xc1,0xc3,0xc6,0xb3,0xb5,0xbc,0xe6,0xf0,0xfa,0xf5,0xf9,0xfd,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xf2,0xf4,0xf7,0xf4,0xf7,0xf8,0xeb,0xee,0xf1,0xda,0xe1,0xe7,0xcc,0xd5,0xdd,0xc8,0xd1,0xd9,0xce,0xd6,0xde,0xd8,0xdd,0xe3,0xdc,0xdf,0xe4,0xec,0xed,0xee,0xee,0xef,0xf0,0xed,0xed,0xee,0xe1,0xe2,0xe4,0xcf,0xd3,0xd7,0xc9,0xcc,0xd2,0xbf,0xc3,0xca,0xb3,0xb8,0xc2,0x7d,0x89,0x9d,0x4d,0x65,0x83,0x49,0x60,0x7f,0x46,0x5d,0x7d,0x47,0x5d,0x7c,0x4a,0x60,0x7d,0x6b,0x79,0x8d,0xc5,0xc6,0xc9,0xaf,0xb2,0xb7,0xf7,0xfa,0xfd,0xf2,0xf6,0xfb,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xf2,0xf4,0xf6,0xd6,0xdd,0xe6,0xdb,0xe1,0xe7,0xe0,0xe5,0xe9,0xef,0xef,0xf1,0xd1,0xd6,0xdc,0xe0,0xe2,0xe6,0xe7,0xe9,0xeb,0xda,0xdd,0xe0,0xd6,0xda,0xdd,0xa4,0xae,0xbc,0x82,0x91,0xa6,0x79,0x87,0x9e,0x6f,0x80,0x98,0x5b,0x70,0x8b,0x52,0x69,0x85,0x4b,0x63,0x81,0x49,0x60,0x7f,0x49,0x60,0x7f,0x48,0x5e,0x7d,0x4e,0x63,0x7f,0xaf,0xb3,0xb9,0x95,0x9a,0xa6,0xfe,0xfe,0xff,0xfe,0xfe,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xf7,0xf8,0xf8,0xe3,0xe8,0xed,0xf1,0xf2,0xf5,0xfa,0xfa,0xfa,0xed,0xee,0xf0,0xdd,0xe0,0xe5,0xeb,0xec,0xee,0xe7,0xe8,0xea,0xe0,0xe3,0xe5,0xde,0xdf,0xe3,0x9d,0xa8,0xb7,0x65,0x7a,0x99,0x5b,0x74,0x91,0x5a,0x72,0x8e,0x54,0x6c,0x89,0x51,0x68,0x86,0x4d,0x65,0x83,0x4b,0x62,0x81,0x4a,0x62,0x80,0x49,0x5e,0x7e,0x49,0x5f,0x7c,0x67,0x74,0x8a,0x67,0x73,0x88,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xf9,0xf9,0xf9,0xee,0xf0,0xf3,0xf9,0xf9,0xf9,0xfd,0xfd,0xfd,0xef,0xf1,0xf3,0xd5,0xda,0xe0,0xeb,0xed,0xef,0xe3,0xe5,0xe8,0xcf,0xd4,0xda,0xbd,0xc3,0xcc,0x83,0x94,0xaa,0x6f,0x84,0x9f,0x62,0x79,0x95,0x5f,0x75,0x92,0x5d,0x74,0x8f,0x56,0x6d,0x8a,0x50,0x68,0x86,0x51,0x68,0x84,0x4d,0x64,0x82,0x4b,0x61,0x7f,0x48,0x5f,0x7d,0x48,0x5d,0x7b,0x4a,0x5e,0x7a,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xfd,0xfd,0xfd,0xfa,0xfa,0xfb,0xff,0xff,0xfe,0xfd,0xfd,0xfd,0xf8,0xf9,0xf9,0xf6,0xf6,0xf6,0xee,0xf1,0xf2,0xda,0xde,0xe3,0xc9,0xcf,0xd7,0xbe,0xc6,0xd1,0x98,0xa6,0xb8,0x8b,0x9c,0xb1,0x72,0x88,0xa2,0x67,0x7d,0x9a,0x64,0x79,0x95,0x5f,0x76,0x91,0x59,0x70,0x8d,0x55,0x6c,0x88,0x53,0x6a,0x86,0x51,0x67,0x83,0x4d,0x64,0x81,0x4b,0x61,0x7e,0x48,0x5e,0x7c,0x49,0x5d,0x7a,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xfe,0xfe,0xfe,0xfd,0xfd,0xfd,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xfe,0xfe,0xfe,0xff,0xff,0xff,0xfc,0xfd,0xfd,0xff,0xff,0xff,0xfe,0xfe,0xfe,0xfb,0xfb,0xfb,0xf6,0xf6,0xf6,0xe2,0xe6,0xe9,0xb4,0xc1,0xcf,0xb1,0xbb,0xc9,0xa7,0xb1,0xc0,0x9a,0xa8,0xb9,0x7d,0x91,0xa9,0x68,0x80,0x9c,0x66,0x7c,0x97,0x61,0x78,0x91,0x5a,0x71,0x8d,0x57,0x6f,0x8a,0x54,0x6a,0x87,0x51,0x68,0x84,0x51,0x66,0x82,0x4c,0x63,0x7f,0x46,0x5d,0x7b,0x53,0x65,0x7e,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xfe,0xfe,0xfe,0xff,0xff,0xff,0xfe,0xfe,0xfe,0xff,0xff,0xff,0xfc,0xfb,0xfc,0xf5,0xf7,0xf9,0xea,0xed,0xf2,0xe7,0xec,0xf1,0xef,0xf2,0xf5,0xfc,0xfc,0xfb,0xfc,0xfc,0xfc,0xee,0xf0,0xf2,0xd0,0xd7,0xde,0xae,0xbb,0xc9,0x83,0x98,0xb0,0x7f,0x94,0xac,0x7c,0x90,0xa9,0x6a,0x82,0x9e,0x68,0x7e,0x9a,0x61,0x78,0x93,0x5d,0x75,0x90,0x59,0x70,0x8b,0x55,0x6c,0x87,0x53,0x69,0x85,0x52,0x68,0x83,0x4f,0x65,0x80,0x48,0x5f,0x7a,0x4b,0x5f,0x79,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xfa,0xfc,0xfe,0xf3,0xf8,0xfd,0xe4,0xf0,0xfa,0xd8,0xe8,0xf7,0xf2,0xf6,0xfa,0xf6,0xf7,0xf9,0xd1,0xdf,0xee,0xd5,0xe0,0xec,0xdf,0xe7,0xef,0xee,0xf1,0xf3,0xf4,0xf5,0xf6,0xe2,0xe6,0xec,0xf1,0xf1,0xf3,0xd1,0xd9,0xe0,0xaa,0xb9,0xca,0xb6,0xc2,0xce,0x82,0x99,0xb0,0x7d,0x94,0xaa,0x6f,0x87,0xa1,0x6e,0x84,0x9e,0x65,0x7c,0x96,0x5e,0x75,0x90,0x5b,0x72,0x8c,0x56,0x6d,0x89,0x54,0x6a,0x86,0x51,0x67,0x82,0x4c,0x62,0x7f,0x47,0x5d,0x7a,0x46,0x5c,0x77,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xfc,0xfd,0xfe,0xf9,0xfb,0xfd,0xd4,0xe5,0xf6,0xcc,0xde,0xf3,0xd0,0xe0,0xf2,0xf0,0xf4,0xf8,0xc9,0xda,0xeb,0xd6,0xe2,0xee,0xfe,0xfe,0xfe,0xee,0xf0,0xf4,0xd8,0xe0,0xe7,0xe9,0xec,0xef,0xd6,0xdc,0xe3,0xe1,0xe5,0xe9,0xd3,0xd9,0xde,0xcc,0xd3,0xda,0xae,0xbb,0xc8,0x9a,0xa8,0xba,0x6f,0x86,0xa0,0x7d,0x90,0xa6,0x7d,0x8e,0xa3,0x5c,0x74,0x90,0x58,0x71,0x8c,0x55,0x6c,0x88,0x52,0x68,0x84,0x4d,0x64,0x80,0x4a,0x60,0x7d,0x47,0x5d,0x78,0x46,0x5b,0x75,0xf0,0xf5,0xfc,0xf2,0xf6,0xfc,0xf8,0xfa,0xfd,0xfa,0xfc,0xfd,0xcf,0xe0,0xf2,0xc4,0xda,0xf0,0xba,0xd2,0xeb,0xc1,0xd5,0xeb,0xc8,0xda,0xec,0xd7,0xe3,0xef,0xf1,0xf4,0xf7,0xf4,0xf5,0xf6,0xbd,0xcf,0xe0,0xaf,0xc4,0xd9,0xde,0xe4,0xe9,0xd7,0xdd,0xe3,0xea,0xeb,0xed,0xcf,0xd5,0xdc,0xc0,0xc8,0xd2,0xbd,0xc5,0xcf,0x7c,0x91,0xaa,0x6c,0x84,0x9f,0x68,0x7d,0x99,0x66,0x7b,0x94,0x5b,0x72,0x8d,0x58,0x6f,0x89,0x52,0x6a,0x84,0x4e,0x65,0x81,0x49,0x61,0x7d,0x4a,0x60,0x7b,0x48,0x5d,0x78,0x46,0x5a,0x74,};
