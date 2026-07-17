/*
  This sketch is written for a Arduino Mega 2560 Rev3 and it is moving a stepper motor with a DM542 microstep driver attached to a TIRF rail. If a digital signal 
  is send from a NI DAQ 6713 card (https://www.ni.com/de-de/support/model.pci-6713.html?srsltid=AfmBOoo-oGrrvV_yyremIEEHGdVB5U6n_pt5i-xaASUSRoXrJi6zpJkq), 
  the stepper motor moves the TIRF rail from TIRF to Epi and back to TIRF configuration. Furthermore a 405 nm diode laser 
  (https://www.toptica.com/products/single-mode-diode-lasers/ibeam-smart) is pulsed if a digital signal is send by the DAQ card.
*/

///////////////
//Define pins//
///////////////

// Define UV diode laser pins:
int UV_laser_digi_pin = 50; // connected to digital input of 405 nm iBeam smart laser
int UV_laser_manual_pin = 51; // connect to digital signal from DAQ card to manual switch on UV laser
int UV640 = 52; // digital signal from DAQ card to activate UV640 pulse program
int UV532 = 53; // digital signal from DAQ card to activate UV532 pulse program

// Define stepper motor pins:
int pul_pin = 5;
int dir_pin = 6;
int enbl_pin = 7;
int trigger_pin = 49;


/////////////////////////////
//define initial parameters//
/////////////////////////////

// phases for pulsing UV laser after red/UV640 program:
unsigned long first_phase = 120000L;
unsigned long second_phase = 120000L;
unsigned long third_phase = 60000L;

// phases for pulsing UV laser after green/UV532 program:
unsigned long green_first_phase = 120000L;
unsigned long green_second_phase = 30000L;

// UV laser:
int on_time = 100;
int off_time = 14900;
unsigned long c_time = millis();
int UV640_state;
int UV532_state;
int UV_manual;

// stepper motor velocity v in microseconds between steps:
int v = 1300; // low expressed 1300; high expressed 5000
// number of steps in one direction:
int  steps = 3500;
int Index;
int previousState = 0;




void setup(){
  Serial.begin(9600);   // open connection to PC for debugging

  // initialize UV-laser:
  pinMode(UV_laser_digi_pin, OUTPUT);
  pinMode(UV_laser_manual_pin, INPUT);
  digitalWrite(UV_laser_digi_pin, LOW);
  pinMode(UV640, INPUT);
  pinMode(UV532, INPUT);

  // initialize stepper motor:
  pinMode(pul_pin, OUTPUT);
  pinMode(dir_pin, OUTPUT);
  pinMode(enbl_pin, OUTPUT);
  pinMode(trigger_pin, INPUT_PULLUP);
  digitalWrite(enbl_pin, HIGH); // deactivate stepper motor --> TIRF rail can be moved manual
}




int move_stepper() {
  // move stepper motor from TIRF to Epi and back to TIRF configuration
  digitalWrite(enbl_pin, LOW); // enable motor control
  delayMicroseconds(6);
  digitalWrite(dir_pin, LOW); // set direction to move from TIRF to Epi
  delayMicroseconds(6);

  for (Index = 0; Index < steps; Index++) // move rail until Epi is reached
  {
    digitalWrite(pul_pin, HIGH);
    delayMicroseconds(v);
    digitalWrite(pul_pin, LOW);
    delayMicroseconds(v);
  }

  delay(500);
  digitalWrite(dir_pin, HIGH); // set direction to move from Epi to TIRF
  delayMicroseconds(6); 

  for (Index = 0; Index < steps; Index++) // move rail until TIRF is reached again
  {
    digitalWrite(pul_pin, HIGH);
    delayMicroseconds(v);
    digitalWrite(pul_pin, LOW);
    delayMicroseconds(v);
  }

  delay(5);
  digitalWrite(enbl_pin, HIGH); // deactivate stepper motor --> TIRF rail can be moved manual
}




void loop(){

  ////////////////////
  // move TIRF rail //
  ////////////////////

  int triggerState = digitalRead(trigger_pin);
  // if stepper motor is triggered wait and check again if it is pulse or steady as on program start
  if (triggerState == HIGH && previousState == LOW){
    // wait until pulse should have finished
    delay(150);

    // check triggerState second time 
    int second_trigger = digitalRead(trigger_pin);
    // pulse -> move or steady as on program start/ not connected -> wait
    if (second_trigger == LOW){
      move_stepper();
    }
    else{
      delay(50);
    }
  }
  previousState = triggerState;


  ////////////////////////////////
  // control and pulse UV laser //
  ////////////////////////////////

  // check if UV laser was manual activated:
  UV_manual = digitalRead(UV_laser_manual_pin);
  if (UV_manual){
    digitalWrite(UV_laser_digi_pin, HIGH);
  }
  else {
    digitalWrite(UV_laser_digi_pin, LOW);
  }


  // check if UV puls program was activated:
  UV640_state = digitalRead(UV640);
  UV532_state = digitalRead(UV532);

  if (UV640_state){
    // 1.Phase: Wait for first phase
    c_time = millis();
    while ((millis() < (first_phase + c_time)) && UV640_state) {
      digitalWrite(UV_laser_digi_pin, LOW);
      UV640_state = digitalRead(UV640);
      Serial.println("1.phase");
    }

    // 2.Phase: Pulse extra low
    c_time = millis();
    while ((millis() < (second_phase + c_time)) && UV640_state) {
      digitalWrite(UV_laser_digi_pin, HIGH);
      delayMicroseconds(100);
      digitalWrite(UV_laser_digi_pin, LOW);
      delay(500);
      UV640_state = digitalRead(UV640);
      Serial.println("2.phase");
    }

    // 3.Phase: Pulse low
    c_time = millis();
    while ((millis() < (third_phase + c_time)) && UV640_state) {
      digitalWrite(UV_laser_digi_pin, HIGH);
      delayMicroseconds(100);
      digitalWrite(UV_laser_digi_pin, LOW);
      delayMicroseconds(14900);
      UV640_state = digitalRead(UV640);
      Serial.println("3.phase");
    }

    // 4.Phase: Ramp
    c_time = millis();
    on_time = 100;
    off_time = 14900;
    while ((on_time < 14900) && UV640_state) {
      digitalWrite(UV_laser_digi_pin, HIGH);
      delayMicroseconds(on_time);
      digitalWrite(UV_laser_digi_pin, LOW);
      delayMicroseconds(off_time);
      on_time += 1;
      off_time -= 1;

      UV640_state = digitalRead(UV640);
      Serial.println("4.phase");
    }

    // 5.Phase: continous
    while (UV640_state){
      digitalWrite(UV_laser_digi_pin, HIGH);
      UV640_state = digitalRead(UV640);
      Serial.println("5.phase");
    }
  }

  if (UV532_state){
    // 1.Phase: Wait for first phase
    c_time = millis();
    while ((millis() < (green_first_phase + c_time)) && UV532_state) {
      digitalWrite(UV_laser_digi_pin, LOW);
      UV532_state = digitalRead(UV532);
      Serial.println("green 1.phase");
    }

    // 2.Phase: Pulse low
    c_time = millis();
    while ((millis() < (green_second_phase + c_time)) && UV532_state) {
      digitalWrite(UV_laser_digi_pin, HIGH);
      delayMicroseconds(100);
      digitalWrite(UV_laser_digi_pin, LOW);
      delayMicroseconds(14900);
      UV532_state = digitalRead(UV532);
      Serial.println("green 2.phase");
    }

    // 3.Phase: Ramp
    c_time = millis();
    on_time = 100;
    off_time = 14900;
    while ((on_time < 14900) && UV532_state) {
      digitalWrite(UV_laser_digi_pin, HIGH);
      delayMicroseconds(on_time);
      digitalWrite(UV_laser_digi_pin, LOW);
      delayMicroseconds(off_time);
      on_time += 1;
      off_time -= 1;

      UV532_state = digitalRead(UV532);
      Serial.println("green 3.phase");
    }

    // 4.Phase: continous
    c_time = millis();
    while (UV532_state) {
      digitalWrite(UV_laser_digi_pin, HIGH);
      UV532_state = digitalRead(UV532);
      Serial.println("green 4.phase");
    }

  }

  digitalWrite(UV_laser_digi_pin, LOW);



  // activate to send signals for debugging:
  //Serial.println("on time");
  //Serial.println(on_time);
  //Serial.println(millis()+wait_time);
  //Serial.println(off_time);
  //Serial.println("");

}
