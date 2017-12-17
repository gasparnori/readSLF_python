import load_intan_rhd_format as intan
import time
import numpy as np
import math
import pprint
import matplotlib.pyplot as plt



#supply voltage is only sampled once per data block
filePath="D:/Nora/amir_recording"
WriteInfos=False
MainsFreqChannel=5
SyncChannel=4
MaxSize=100000       #maximum sample size to be read from the files

class SLFData:
    def __init__(self, folder):
        self.folderName=folder
        self.readInfo()
        self.readtimeData()
        self.readDigital()
        self.readAllInOne()
        self.SupplyVoltage()

    def readInfo(self):
        res=intan.read_data(self.folderName+'/info.rhd')
        self.info=IntanInfos(res)

    def readtimeData(self):
        #data was sampled 60 times per data block
        self.timedata = np.fromfile(self.folderName+'/time.dat',dtype=np.int32) #int32!!!
        self.num_samples = len(self.timedata)
        self.dt=1/self.info.amplifier_sampling_rate
        self.scaledTimes=self.timedata*self.dt

        #because the supply voltage was only sampled once in a data block
        self.supplytime=self.timedata[::60]
        #auxiliary inputs were only sampled 15 times per block          -------------not sure if needed
        self.auxtime=self.timedata[::4]
        if WriteInfos:
            print ('Total recording time: ', self.scaledTimes[-1]-self.scaledTimes[0])
            print ('number of samples', self.num_samples)

    def readAllInOne(self):
        self.AllInOnedata = np.fromfile(self.folderName+'/all_in_one.dat',dtype=np.int16) #int16!!!!
        self.numAllInSamples=len(self.AllInOnedata)
        sumChannels=self.info.NumAuxCh+self.info.NumADCCh+self.info.NumAmpCh


        #####read channels from the all in one data with parameters: offset, number of channels, number of all channels
        self.ampData=self.readChannelData(0, self.info.NumAmpCh, sumChannels)
        self.ADCData=self.readChannelData((self.info.NumADCCh), self.info.NumADCCh, sumChannels,)
        self.auxData= self.readChannelData((self.info.NumAmpCh+self.info.NumADCCh), self.info.NumAuxCh, sumChannels)
        if WriteInfos:
            print ('all in one data read: \n   number of smplifier channels: ', self.info.NumAmpCh, "\n   number of ADC channels: ", self.info.NumADCCh, "\n   number of auxiliary channels: ", self.info.NumAuxCh)
        ############################### Scale voltage levels
        ############# see Intan RHD2000 datasheet\electrical characteristics
        if self.ampData is not None:
            for i in range(0, len(self.ampData)):
                self.ampData[i] = [ 0.195 * float(k) for k in self.ampData[i]] #[uV]
        if self.ADCData is not None:
            for i in range(0, len(self.ADCData)):
                self.ADCData[i] = [ (10/pow(2,16)) * float(k) for k in self.ADCData[i]] #[-5V to 5V].
                #print(self.ADCData)
        if self.auxData is not None:
            for i in range(0, len(self.auxData)):
                self.auxData = [(2.45/2^16) * (float(k) + 2^15) for k in self.auxData[i]] # [0.1V to 2.45 V];

    def readChannelData(self, offset, num_channels, sum_channels):
        if num_channels>0:
            d=[]
            for i in range(0, num_channels):
                minRange=offset+i
                #print(minRange)
                if MaxSize>((len(self.AllInOnedata)-offset)/sum_channels):
                    maxRange=len(self.AllInOnedata)
                else:
                    maxRange=(MaxSize+offset+i)*sum_channels
                #print(maxRange)
                d.append(self.AllInOnedata[minRange:maxRange:sum_channels])
                #print(d)
            return d
        else:
            return None

    def readDigital(self):
        self.digitData=self.DChannel(self.info.NumdigCh, self.folderName)

    def SupplyVoltage(self):
        self.supply=np.fromfile(self.folderName+'/supply.dat',dtype=np.uint16) #int16!!!
        if self.supply is not None:
            self.supply = 74.8e-6 * self.supply # [V]
            if WriteInfos:
                print ("supply voltage read too", len(self.supply))

    def DChannel(self, numChannels, folderName):
        rawdata = np.fromfile(folderName+'/digitalin.dat',dtype=np.uint16) #int16!!!

    # Digital data is saved as 16-bit words in uint16 format. For example,
    # if digital inputs 0, 4, and 5 are high and the rest low, the uint16
    # value for this sample time will be 2^0 + 2^4 + 2^5 = 1 + 16 + 32 =
    # 49. The following is an adaptation of MATLAB XXX function. From right
    # to left: raw data (digitalin) is converted to a matrix with each
    # channel (column) devided by 2^channel number. The values are rounded
    # to the nearest low integer and converted to binary by calculating the
    # reminder with 2 (odd - 1; even - 0). Alas, the matrix is flipped such
    # that the least significant number (ch 1) is on the right. For
    # comprehension, run segments of the following line with the value 8209
    # (which represents channels 1, 4 and 13 high, the rest low) instead of
    # digitalin.
        if rawdata is not None:
            data=[]
            for i in range (0,15):
                d=np.zeros(MaxSize)
                for j in range (0, MaxSize):
                    d[j]=(rawdata[j] >> (i)) & 1
                data.append(d)
            return data
        else:
            return None

        ################################Nora's comment: still not sure if it works perfectly, use with care!!!!
    def PlotAllDigital(self):
        plt.title('All Digital channels')
        for i in range(0,self.info.NumdigCh):
            plt.subplot(math.ceil(self.info.NumdigCh/3),3,i+1)
            plt.plot(self.scaledTimes[0:MaxSize], self.digitData[i][0:MaxSize])
            plt.ylabel("channel "+str(i))
            plt.xlabel("scaled times in s")
            plt.axis([0, self.scaledTimes[MaxSize], 0,1])
        plt.show()

    def PlotAllADC(self):
        plt.title('All ADC channels')
        for i in range(0,self.info.NumADCCh):
            plt.subplot(math.ceil(self.info.NumADCCh/3),3,i+1)
            plt.plot(self.scaledTimes[0:MaxSize], self.ADCData[i][0:MaxSize])
            plt.ylabel("channel "+str(i))
            plt.xlabel("scaled times in s")
            plt.axis([0, self.scaledTimes[MaxSize], -5,5])
        plt.show()

    def PlotAllAmplifiers(self):
        plt.title('All Amplifier channels')
        for i in range(0,self.info.NumAmpCh):
            plt.subplot(math.ceil(self.info.NumAmpCh/3),3,i+1)
            plt.plot(self.scaledTimes[0:MaxSize], self.ampData[i][0:MaxSize])
            plt.ylabel("channel "+str(i))
            plt.xlabel("scaled times in s")
            #plt.axis([0, self.scaledTimes[MaxSize], -5,5])
        plt.show()

    def PlotAmplifier(self, channel):
        plt.plot(self.scaledTimes[0:MaxSize], self.ampData[channel])
        plt.ylabel("Amplifier channel "+str(channel))
        plt.xlabel("scaled times in s")
        #plt.axis([0, MaxSize, 0, max(self.digitData[0])])
        plt.show()

    def PlotADC(self, channel):
        plt.plot(self.scaledTimes[0:MaxSize], self.ADCData[channel][0:MaxSize])
        plt.ylabel("ADC Data")
        plt.xlabel("scaled times in s")
        plt.axis([0, self.scaledTimes[MaxSize], -5,5])
        plt.show()

    def PlotDigital(self, channel):
        plt.plot(self.scaledTimes[0:MaxSize], self.digitData[channel])
        plt.ylabel("digital Data")
        plt.xlabel("scaled times in s")
        #plt.axis([0, MaxSize, 0, max(self.digitData[0])])
        plt.show()

    def RemoveLineSignal(self): #############################################still missing
        return

class IntanInfos:
    def __init__(self, infos):
        if WriteInfos:
            pprint.pprint(infos)
        self.notes=infos['notes']
        self.freq_params=infos['frequency_parameters']
        self.amplifier_sampling_rate=self.freq_params['amplifier_sample_rate']
        self.adc_sampling_rate=self.freq_params['board_adc_sample_rate']

        if 'amplifier_channels' in infos:
            self.NumAmpCh =len(infos['amplifier_channels'])
        else:
            self.NumAmpCh=0
        if 'board_adc_channels' in infos:
            self.NumADCCh = len(infos['board_adc_channels'])
        else:
            self.NumADCCh=0
        if 'board_dig_in_channels'  in infos:
            self.NumdigCh = len(infos['board_dig_in_channels'])
        else:
            self.NumdigCh=0
        if 'aux_input_channels' in infos:
            self.NumAuxCh=len(infos['aux_input_channels'])
        else:
            self.NumAuxCh=0
        if 'supply_voltage_channels' in infos:
            self.NumSupplyCh=len(infos['supply_voltage_channels'])
        else:
            self.NumSupplyCh=0



print("Reading Stark Lab format (SLF) Data files")
tic= time.time()
d=SLFData(filePath)
#for i in range(0, len(d.ADCData)):
#    print (d.ADCData[i])
#    d.ADCData[i] = [ (10/pow(2,16)) * float(k) for k in d.ADCData[i]] #[-5V to 5V].
#    print ("\n\n\n\n\n after")
#    print (d.ADCData[i])
#d.PlotDigital(5)
#for i in range(0, len(d.ADCData)):
d.PlotAllDigital()
d.PlotAllADC()
d.PlotAllAmplifiers()
#d.PlotAmplifier(0)
#print (d.scaledTimes)








# def readChannelData(datain, offset, num_channels, sum_channels):
#         if num_channels>0:
#             d=[]
#             for i in range(0, num_channels):
#                 minRange=offset+i
#                 #print(minRange)
#                 if MaxSize>((len(datain)-offset)/sum_channels):
#                     maxRange=len(datain)
#                 else:
#                     maxRange=(MaxSize+offset+i)*sum_channels
#                 #print(maxRange)
#                 d.append(datain[minRange:maxRange:sum_channels])
#                 #print(d)
#             return d
#         else:
#             return None

#data=[1,2,3,4,5,6,7,8,9,10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]
#print(readChannelData(data, 2, 2, 2))

#AllInOnedata = np.fromfile('D:/Nora/amir_recording/all_in_one.dat',dtype=np.int16) #int16!!!!

