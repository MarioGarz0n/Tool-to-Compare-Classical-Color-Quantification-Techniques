#!/bin/bash

RUTA="imagenes/" #ruta para llegar a la carpeta de las imagenes

if [[ $# < 1 ]]
then
   echo "Introduzca la opción que quiera."
   echo "1->Binary-splitting"
   echo "2->K-Means"
   echo "3->Median"
   echo "4->Wu"
   echo "5->Octree"
   echo "6->Variance-Based"
   echo "7->Neuquant"
   exit 1
fi

if [[ $1 = 1 ]]
then
   make bs
   ALGO="BS"
   PROG="./ejec_bs"
elif [[ $1 = 2 ]]
then
   make km
   ALGO="KM"
   PROG="./ejec_kmeans"
elif [[ $1 = 3 ]]
then
   make median
   ALGO="MC"
   PROG="./ejec_median"
elif [[ $1 = 4 ]]
then
   make wu
   ALGO="WU"
   PROG="./ejec_wu"
elif [[ $1 = 5 ]]
then
   make octree
   ALGO="OC"
   PROG="./ejec_octree"
elif [[ $1 = 6 ]]
then
   make vb
   ALGO="VB"
   PROG="./ejec_vb"
elif [[ $1 = 7 ]]
then
   make neuquant
   ALGO="NQ"
   PROG="./ejec_nq"
else
   echo "Introduzca la opción que quiera."
   echo "1->Binary-splitting"
   echo "2->K-Means"
   echo "3->Median"
   echo "4->Wu"
   echo "5->Octree"
   echo "6->Variance-Based"
   echo "7->Neuquant"
   exit 2
fi

clear

# Rótulo para los 4 ficheros que almacenarán los errores para cada tamaño de paleta
for C in 32 64 128 256 ;
do 
   echo "N Img MSE PSNR MAE BRISQUE DSS FSIM GMSD HaarPSI MDSI MS-GMSDc SSIM ms-SSIM iw-SSIM TV-index VIFP VSI SR-SIM TIEMPO" > IQI_${ALGO}_${C}.txt
   echo "N Img MSE_${ALGO} PSNR_${ALGO} MAE_${ALGO} BRISQUE_${ALGO} DSS_${ALGO} FSIM_${ALGO} GMSD_${ALGO} HaarPSI_${ALGO} MDSI_${ALGO} MS-GMSDc_${ALGO} SSIM_${ALGO} ms-SSIM_${ALGO} iw-SSIM_${ALGO} TV-index_${ALGO} VIFP_${ALGO} VSI_${ALGO} SR-SIM_${ALGO} TIEMPO_${ALGO}" >> IQI_${ALGO}_${C}.txt
done


# Para cada tamaño de paleta cuantizada
for C in 32 64 128 256 ;
do	
	#para cada imagen
	for F in adirondack_chairs.ppm astro_bodies.ppm astronaut.ppm balinese_dancer.ppm ball_caps.ppm birthday_baloons.ppm bosnian_pine_needle.ppm buggy.ppm calaveras.ppm carrots.ppm chalk_pastels.ppm chicken_dish.ppm chili_peppers.ppm clownfish.ppm color_chart.ppm  color_checker.ppm coloring_pencils.ppm columbia_crew.ppm common_jezebel.ppm common_lantanas.ppm  cosmic_vista.ppm  craft_cards.ppm crepe_paper.ppm cruise_ship.ppm curler.ppm daisy_bouquet.ppm daisy_poster.ppm  easter_egg_basket.ppm  easter_eggs.ppm eastern_rosella.ppm felt_ball_trivet.ppm fishing_nets.ppm floating_market.ppm fruit_dessert.ppm fruit_stand.ppm  fruits.ppm fruits_veggies.ppm german_hot_air_balloon.ppm girl.ppm gourds.ppm grilled_food.ppm  hard_candy.ppm italian_hot_air_balloon.ppm jacksons_chameleon.ppm king_penguin.ppm king_vulture.ppm  kingfisher.ppm korean_dancer.ppm lights.ppm macarons.ppm macaws.ppm malayan_banded_pitta.ppm  mandarin_ducks.ppm mandarinfish.ppm mangoes.ppm marrakech_museum.ppm maya_beach.ppm  medicine_packets.ppm moroccan_babouches.ppm motocross.ppm motorcycle.ppm mural.ppm  nylon_cords.ppm paper_clips.ppm peacock.ppm pencils.ppm pigments.ppm  pink_mosque.ppm plushies.ppm prickly_pears.ppm puffin.ppm race_car.ppm  red_eyed_tree_frog.ppm red_knobbed_starfish.ppm rescue_helicopter.ppm rose_bouquet.ppm sagami_temple.ppm salad_bowl.ppm schoolgirls.ppm seattle_great_wheel.ppm shawls.ppm shopping_bags.ppm  siberian_tiger.ppm skiers.ppm spices.ppm sports_bicycles.ppm sun_parakeet.ppm tablet.ppm  textile_market.ppm trade_fair_tower.ppm traffic.ppm tulip_field.ppm umbrellas.ppm  veggie_pizza.ppm veggies.ppm venetian_lagoon.ppm vintage_cars.ppm wooden_toys.ppm  wool_carder_bee.ppm yasaka_pagoda.ppm ;
	do
		# NÚMERO DE TEST INDEPENDIENTES EJECUTADOS PARA UNA MISMA CONFIGURACIÓN
 		# -- NECESARIO SOLO PARA MÉTODOS NO DETERMINÍSTICOS (es decir, los que no generan la misma imagen con los mismos parámetros) --
		for TEST in {1..20};
		do 
			#ejecuto el programa, con el tamaño de la paleta indicado y se vuelca la salida 
			$PROG ${RUTA}${F} ${C} >> dump_${ALGO}.txt
			
			#Calculo múltiples medidas de error sobre la imagen cuantizada que acabo de generar y vuelco los valores a un fichero
			python3 metricas_error.py ${RUTA}${F} ${ALGO}_${C}_${F} ${C} >> IQI_${ALGO}_${C}.txt
			
			#Borro la imagen nueva, que ya no necesito (guardarlas todas ocupa mucho disco duro)
			rm ${ALGO}_${C}_${F}
			
			#si no es el kmeans se tendrán los mismos resultados todas las veces (deterministicos)
			if [[ $1 != 2 ]]
			then
				break
			fi			
		done #tests sucesivos
	done #imagen original
done #colores de la paleta cuantizada

#se calcula la media del tiempo para el kmeans
if [[ $1 = 2 ]] 
then
	bash media_tiempo.sh dump_${ALGO}.txt
	cat auxiliar_media_tiempo.txt > dump_${ALGO}.txt
	rm auxiliar_media_tiempo.txt
fi	

#se guardan los errores y los tiempos en libros de cálculo
cat dump_${ALGO}.txt > tiempo_${ALGO}.ods
rm dump_${ALGO}.txt 

for C in 32 64 128 256 ; 
do
	if [[ $1 = 2 ]] #media de errores para kmeans, el unico no deterministico
	then
		bash media_errores.sh IQI_${ALGO}_${C}.txt
		cat archivo_auxiliar_errores.txt > IQI_${ALGO}_${C}.txt
	   rm archivo_auxiliar_errores.txt
	fi
	
	cat IQI_${ALGO}_${C}.txt > errores_${ALGO}_${C}.ods
	rm IQI_${ALGO}_${C}.txt
done








