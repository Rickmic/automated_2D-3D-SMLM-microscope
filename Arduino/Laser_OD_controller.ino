/*
  This sketch is written for a Arduino Mega 2560 Rev3 and it modulates a gem 640 and gem 532 DPSS laser
  (https://novantaphotonics.com/product/gem-solid-state-continuous-wave-lasers/) based on the analogue signal from a NI DAQ 6713 card 
  (https://www.ni.com/de-de/support/model.pci-6713.html?srsltid=AfmBOoo-oGrrvV_yyremIEEHGdVB5U6n_pt5i-xaASUSRoXrJi6zpJkq):
  640 nm laser = red_laser
  532 nm laser = green_laser
  GEM640 connected to Serial1 (RX:19, TX:18)
  GEM532 connected to Serial2 (RX:17, TX:16)

  Additionally a motorized filter wheel from Thorlabs (https://www.thorlabs.com/thorproduct.cfm?partnumber=FW103/M) with a BSC201 controller 
  (https://www.thorlabs.com/thorproduct.cfm?partnumber=BSC201) is moved according to the digital signal from the DAQ card. Furthermore 
  a servo motor with an attached OD-wheel is moved according to a analogue signal from a NI DAQ 6713 card.
*/
#include <SoftwareSerial.h>
#include <Servo.h>

// activate Servo
Servo OD_Servo;

///////////////
//Define pins//
///////////////

// Define DPSS laser pins:
int red_laser_pin = A15; // analogue signal from DAQ card to modulate red laser
int green_laser_pin = A14; // analogue signal from DAQ card to modulate green laser

// Define filter wheel pins (0: activated, 1: deactivated):
int home_pin = 22; // connected to Pin 4 (Trigger In) on BSC201, Hardware Triggering of BSC201 must be activated and set to "Home" as action.
int fwd_pin = 23; // connected to Pin 2 (Jog forwards) on BSC201
int bwd_pin = 24; // connected to Pin 10 (Jog backwards) on BSC201
int filterW_pin_A = A12; // analogue signal A from DAQ card to control filter wheel
int filterW_pin_B = A13; // analogue signal B from DAQ card to control filter wheel

// Define Servo pins:
int servo_pin = 8;
int OD_wheel_slider = A1;
int new_OD_pos;

/////////////////////////////
//define initial parameters//
/////////////////////////////

// set laser state to off:
int red_laser_state = 0;
int green_laser_state = 0;

// set initial laser power to 0 mW:
int red_laser_power = 0;
int green_laser_power = 0;
int new_red_laser_power = 0;
int new_green_laser_power = 0;

// set initial FW position to home:
int FW_pos = 1;

void setup(){
  // initialize laser:
  Serial.begin(9600);   // open connection to PC for debugging
  Serial1.begin(57600); // initialize serial port for GEM640 laser (Baud-Rate 57600)
  Serial2.begin(9600);  // initialize serial port for GEM532 laser (Baud-Rate 9600)
  Serial1.write("\r");  // clean up 
  Serial2.write("\r");
  delay(50);
  Serial1.write("OFF\r");     // GEM640 is switched initially off
  Serial2.write("OFF\r");     // GEM532 is switched initially off
  delay(50);

  // initialize filter wheel:
  pinMode(fwd_pin, OUTPUT);
  pinMode(bwd_pin, OUTPUT);
  pinMode(home_pin, OUTPUT);
  pinMode(filterW_pin_A, INPUT_PULLUP);
  pinMode(filterW_pin_B, INPUT_PULLUP);
  digitalWrite(fwd_pin, 1);
  digitalWrite(bwd_pin, 1);
  delay(100);
  digitalWrite(home_pin, 0);
  delay(300);
  digitalWrite(home_pin, 1);
  delay(5000); // wait until filter wheel is homed

  // initialize OD wheel servo:
  pinMode(OD_wheel_slider, INPUT_PULLUP);
  OD_Servo.attach(servo_pin);
  OD_Servo.writeMicroseconds(1500);
  
  
}

int calcLaserPower(int slider){
  // calcLaserPower calculates the laser power from the analogue slider value
  int laser_power;
  if (slider == 0) { 
    laser_power = 0;
  }
  else {
    laser_power = slider/2.046;
  }
  return laser_power;
}

int moveFWposition(int FW_pos, int new_pos){
  // moveFWposition moves Filter Wheel to a new position depending on the old Filter Wheel position
  int difference;
  difference = new_pos - FW_pos;
  if (difference > 0){ // move reverse
    while (difference != 0){
      digitalWrite(bwd_pin, 0);
      delay(350);
      digitalWrite(bwd_pin, 1);
      delay(350);
      difference -= 1;
    }
  }
  if (difference < 0){ // move foreward
    while (difference != 0){
      digitalWrite(fwd_pin, 0);
      delay(350);
      digitalWrite(fwd_pin, 1);
      delay(350);
      difference += 1;
    }
  }
}


void loop(){
  // read laser power:
  int red_laser_slider = analogRead(red_laser_pin);
  int green_laser_slider = analogRead(green_laser_pin);

  // calculate new laser power:
  new_red_laser_power = calcLaserPower(red_laser_slider);
  new_green_laser_power = calcLaserPower(green_laser_slider);

  ////////////////////////////////
  // communicate with red laser://
  ////////////////////////////////
  if (new_red_laser_power != red_laser_power && new_red_laser_power > 10) { //red laser power was adjusted and is more then 10 mW
    if (red_laser_state == 0) { // red laser is switched off
      Serial1.write("ON\r"); // switch on red laser on before adjusting laser power
      delay(50);
    }
    // Adjust red laser power:
    red_laser_power = new_red_laser_power;
    Serial1.write("POWER=");
    delay(50);
    Serial1.print(red_laser_power);
    delay(50);
    Serial1.write("\r");
    red_laser_state = 1;
  }
  else if (new_red_laser_power != red_laser_power && new_red_laser_power <= 10) { // red laser power was adjusted and is less then 10 mW
    red_laser_power = new_red_laser_power;
    Serial1.write("OFF\r"); 
    red_laser_state = 0;
    delay(50); 
  }

  //////////////////////////////////
  // communicate with green laser://
  //////////////////////////////////
  if (new_green_laser_power != green_laser_power && new_green_laser_power > 10) { //green laser power was adjusted and is more then 10 mW
    if (green_laser_state == 0) { // green laser is switched off
      Serial2.write("ON\r"); // switch on green laser on before adjusting laser power
      delay(50);
    }
    // Adjust laser power:
    green_laser_power = new_green_laser_power;
    Serial2.write("POWER=");
    delay(50);
    Serial2.print(green_laser_power);
    delay(50);
    Serial2.write("\r");
    green_laser_state = 1;
  }
  else if (new_green_laser_power != green_laser_power && new_green_laser_power <= 10) { // green laser power was adjusted and is less then 10 mW
    green_laser_power = new_green_laser_power;
    Serial2.write("OFF\r"); 
    green_laser_state = 0;
    delay(50); 
  }

  ///////////////////////////////////
  // communicate with filter wheel://
  ///////////////////////////////////

  // read filter wheel input:
  int FW_pin_A_state = digitalRead(filterW_pin_A);
  int FW_pin_B_state = digitalRead(filterW_pin_B);
 
  if (not FW_pin_A_state && not FW_pin_B_state){ // move to filter wheel position 1
    //Serial.println("pos_1");
    moveFWposition(FW_pos, 1);
    FW_pos = 1;
  }
  else if (not FW_pin_A_state && FW_pin_B_state){ // move to filter wheel position 1
    //Serial.println("pos_2");
    moveFWposition(FW_pos, 2);
    FW_pos = 2;
  }
  else if (FW_pin_A_state && not FW_pin_B_state){ // move to filter wheel position 1
    //Serial.println("pos_3");
    moveFWposition(FW_pos, 3);
    FW_pos = 3;
  }
  else if (FW_pin_A_state && FW_pin_B_state){ // move to filter wheel position 1
    //Serial.println("pos_4");
    moveFWposition(FW_pos, 4);
    FW_pos = 4;
  }

  /////////////////
  //move OD-wheel//
  /////////////////

  int OD_slider_val = analogRead(OD_wheel_slider); // read Servo Pin from NIDAQ
  // change OD Wheel pos
  new_OD_pos = map(OD_slider_val, 14, 1023, 1500, 1650); // original 1670
  OD_Servo.writeMicroseconds(new_OD_pos);



  // activate to send signals for debugging:
  Serial.println("red laser power [mW]:");
  Serial.println(red_laser_power);
  Serial.println("green laser power [mW]:");
  Serial.println(green_laser_power);

  //Serial.println("pin_A");
  //Serial.println(FW_pin_A_state);
  //Serial.println("pin_B");
  //Serial.println(FW_pin_B_state);

  delay(100);
}
